! $Id: checkpt_mod.f,v 1.23 2012/04/25 22:46:23 nicolas Exp $
      MODULE CHECKPT_MOD
!
!******************************************************************************
!  Module CHECKPT_MOD contains variables and routines which are used to read
!  and write GEOS-CHEM checkpoint files, which contain tracer concentrations
!  in [v/v] mixing ratio, humidities, temperatures and exit values from rpmares
!  (dkh, 8/27/04, adj_group 6/09/09)
!
!  Module Variables:
!  ============================================================================
!  (1 ) INPUT_CHECKPT_FILE   : Full path name of the checkpt file to be read
!  (2 ) OUTPUT_CHECKPT_FILE  : Full path name (w/ tokens!) of output file
!  (3 ) INPUT_OBS_FILE       : Full path name of the obs file to be read
!  (4 ) OUTPUT_OBS_FILE      : Full path hname (w/tokens!) of obs file 
!
!  Module Routines:
!  ============================================================================
!  (1 ) MAKE_CHECKPT_FILE    : Writes checkpoint file to disk 
!  (2 ) READ_CHECKPT_FILE    : Reads checkpoint file from disk 
!  (3 ) READ_OBS_FILE	     : Read obs file from disk  (include this here 
!                              as the observation file is currently the same
!                              as the checkpt file)
!
!  GEOS-CHEM modules referenced by restart_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f : Module containing routines for binary punch file I/O
!  (2 ) error_mod.f : Module containing NaN and other error check routines
!  (3 ) file_mod.f  : Module containing file unit numbers and error checks
!  (4 ) grid_mod.f  : Module containing horizontal grid information
!  (5 ) time_mod.f  : Module containing routines for computing time & date
!  (6 ) restart_mod.f : Module containing CHECK_DIMENSIONS
!
!  NOTES:
!   Pretty much like a stripped down version of RESTART_MOD  (dkh,8/30/04)
!  (2 ) Swtich from OBS and RP_OUT to using OBS_STT and CHK_STT
!  (3 ) Add CHK_PSC.  (dkh, 03/16/05)
!  (4 ) Added support for full chemistry.  Add module varialbe PART_CASE. 
!       Added subroutine CHECK_DIMENSIONS_2.  Modified READ / WRITE CHK 
!       routines and INIT / CLEAN to support full chem.. 
!       Add SMVGARRAY. 
!       (dkh, 07/22/05)
!  (5 ) Add support for sulfate chemistry -- add SO2_CHK and H2O2_CHK.
!       (dkh, 10/12/05)
!       add WETD_CHK_H2O2s_CHEMT, WETD_CHK_H2O2s_DYNT, etc. (dkh, 10/23/05) 
!  (6 ) Add WETD_CHK_SO2_CHEMT and WETD_CHK_SO2_DYNT. (dkh, 10/31/05)  
!  (7 ) Add CONV_CHK_H2O2s_CHEMT and CONV_CHK_SO2s_DYNT. (dkh, 11/22/05)  
!  (8 ) Add routines MAKE_SAVE_FILE and EXPAND_NAME. (dkh, 07/19/06)   
!  (9 ) Add SOILNOX_CHK. (dkh, 02/06/07) 
!  (10) Add CHK_STT_CON(:,:,:) array for checkpointing STT before convection.
!       (mak, 8/2/07)
!  (11) Add MAKE_CHK_DYN_FILE and READ_CHK_DYN_FILE, move checkpointing of 
!        variables that change at dynamic time steps to these routines. (dkh, 02/02/09) 
!  (12) Update to v8, delete obsolete arrays (like CHK_STT) or ones that are somewhere else
!        now (adj_group, 6/09/09)
!  (13) Add support for LADJ_STRAT (hml, dkh, 02/14/12, adj32_025) 
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================    
      CHARACTER(LEN=255) :: INPUT_CHECKPT_FILE  
      CHARACTER(LEN=255) :: OUTPUT_CHECKPT_FILE 
      CHARACTER(LEN=255) :: INPUT_OBS_FILE  
      CHARACTER(LEN=255) :: OUTPUT_OBS_FILE 

      ! Allocatable checkpoint variables
      TYPE (XPLEX), ALLOCATABLE :: RP_IN(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: RP_OUT(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: CHK_STT_CON(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: CHK_STT(:,:,:,:)
      ! move to adj_arrays_mod.f (mak, 6/14/09)
      !TYPE (XPLEX), ALLOCATABLE :: OBS_STT(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: CHK_PSC(:,:,:)
      COMPLEX*16, ALLOCATABLE:: D_CHK_PSC(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: gamaan_fwd(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: gamold_fwd(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: wh2o_fwd(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: ynh4_fwd(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: eror_fwd(:,:,:,:)
      INTEGER, ALLOCATABLE :: exit_fwd(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: gamana_fwd(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: gamas1_fwd(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: gamas2_fwd(:,:,:,:)
      INTEGER, ALLOCATABLE :: nitr_max(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: ORIG_STT(:,:,:,:)
      INTEGER, ALLOCATABLE :: PART_CASE(:)
      TYPE (XPLEX), ALLOCATABLE :: CHK_STT_BEFCHEM(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: CHK_HSAVE(:,:,:)
      COMPLEX*16, ALLOCATABLE :: D_CHK_HSAVE(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: SO2_CHK(:,:,:)
      COMPLEX*16, ALLOCATABLE :: D_SO2_CHK(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: H2O2_CHK(:,:,:)
      COMPLEX*16, ALLOCATABLE :: D_H2O2_CHK(:,:,:) 
      TYPE (XPLEX), ALLOCATABLE :: WETD_CHK_H2O2s(:,:,:)
      COMPLEX*16, ALLOCATABLE :: D_WETD_CHK_H2O2s(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: WETD_CHK_SO2s(:,:,:)
      COMPLEX*16, ALLOCATABLE :: D_WETD_CHK_SO2s(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: WETD_CHK_SO4(:,:,:)
      COMPLEX*16, ALLOCATABLE :: D_WETD_CHK_SO4(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: WETD_CHK_SO2(:,:,:)
      COMPLEX*16, ALLOCATABLE :: D_WETD_CHK_SO2(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: CONV_CHK_H2O2s(:,:,:)
      COMPLEX*16, ALLOCATABLE :: D_CONV_CHK_H2O2s(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: CONV_CHK_SO2s(:,:,:)
      COMPLEX*16, ALLOCATABLE :: D_CONV_CHK_SO2s(:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: SOILNOX_CHK(:,:)
      !TYPE (XPLEX), ALLOCATABLE :: CHK_STT_TD(:,:,:,:)
      TYPE (XPLEX), ALLOCATABLE :: CHK_STT_TC(:,:,:,:)
      !>>>                
      ! Now include adjoint of F (dkh, 10/03/08) 
      TYPE (XPLEX), ALLOCATABLE :: QC_SO2_CHK(:,:,:,:)
      COMPLEX*16, ALLOCATABLE :: D_QC_SO2_CHK(:,:,:,:)
      !<<<

      ! adj_group:  add for checkpointing lightning NOx emissions 
      TYPE (XPLEX),  ALLOCATABLE :: SLBASE_CHK(:,:,:)

      INTEGER, PARAMETER :: NRPIN  = 7 
      INTEGER, PARAMETER :: NRPOUT = 9 
      INTEGER, PARAMETER :: NNNMAX = 50


      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_CHECKPT_FILE( YYYYMMDD, HHMMSS )
!
!******************************************************************************
!  Subroutine MAKE_CHECKPT_FILE creates GEOS-CHEM checkpt files of tracer 
!  mixing ratios (v/v), temp, rh and exit values in binary punch file format. 
!  (dkh, 8/27/04)!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a checkpoint file       
!
!  Passed via CMN:
!  ============================================================================
!  (1 ) TAU    : TAU value (elapsed hours) at start of diagnostic interval
!
!  Passed via ???:
!  ============================================================================
!  (1 ) CHECKPT	: Array of quantities to be checkpointed     
!                               dim=(IIPAR,JJPAR,LLPAR,NCHECKPT) 
!
!  NOTES:
!   Just like MAKE_RESTART_FILE except
!	- only include quantities used as input to RPMARES
!       - include hhmmss in file name
!       - writes files to ADJ_DIR and can zip them
!       dkh, 9/30/04
!  (2 ) Zip *.chk.* files one day at a time in a parallel loop. Add access
!        to GET_TS_CHEM.   (dkh, 11/22/04)  
!  (3 ) Add support for L_RECOMP option to recompute (rather than checkpoint)
!        variables RP_OUT etc.  (dkh, 02/09/05)
!  (4 ) Now write values from CHK_STT. (dkh, 03/03/05)
!  (5 ) Add CHK_PSC. (03/16/05)
!  (6 ) Added support for full chemistry.  Add references to NVAR, CSPEC, IXSAVE,
!       IYSAVE, IZSAVE, NTLOOP_FORKPP, NSRCX
!       Add variables PART_CASE, JLOOP.  Disable L_RECOMP = FALSE option for now. 
!       Added SMVGARRAY. 
!       (dkh, 07/22/05)
!  (7 ) Add SO2_CHK and H2O2_CHK.  (dkh, 10/12/05)
!  (8 ) Add WETD_CHK_H2O2s_CHEMT, WETD_CHK_H2O2s_DYNT, etc. (dkh, 10/23/05) 
!  (9 ) Add WETD_CHK_SO2_CHEMT, WETD_CHK_SO2_DYNT. (dkh, 10/31/05)  
!  (10) Add SOILNOX_CHK. (dkh, 02/06/07) 
!  (11) Now completely split dynamic from chemical time step checkpoints (dkh, 02/01/09) 
!  (12) Remove obsolete options (L_DEL_CHECKPT, L_ZIP_CHECKPT, L_RECOMP), 
!        check for aeroosl simulation (LSULF) and update names to v8 (dkh, 06/11/09) 
!  (13) Now checkpoint XYLAI (dkh, 10/14/09) 
!  (14) BUG FIX: LVARTROP treated correctly (dkh, 01/26/11) 
!  (15) Add support for CH4 simulation (kjw, dkh, 02/12/12, adj32_023) 
!******************************************************************************
!     
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR 
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,          ONLY : EXPAND_DATE, GET_TAU
      USE TRACER_MOD,        ONLY : N_TRACERS            
      USE LOGICAL_MOD,       ONLY : LCHEM , LSULF 
      USE LOGICAL_MOD,       ONLY : LSOILNOX, LLIGHTNOX 
      USE LOGICAL_MOD,       ONLY : LPRT
      USE LOGICAL_ADJ_MOD,   ONLY : LAERO_THERM 
      USE TRACER_MOD,        ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,        ONLY : ITS_AN_AEROSOL_SIM
      USE TRACER_MOD,        ONLY : N_TRACERS
      USE COMODE_MOD,        ONLY : CSPEC_PRIOR , JLOP
      USE GCKPP_ADJ_GLOBAL,  ONLY : NTT
      USE TRACER_MOD,        ONLY : ITS_A_CH4_SIM, STT

     
      ! LVARTROP support for adj (dkh, 01/26/11)
      USE COMODE_MOD,        ONLY : CSPEC_FULL_PRIOR
      USE LOGICAL_MOD,       ONLY : LVARTROP
      USE COMODE_MOD,        ONLY : ISAVE_PRIOR
      USE COMODE_MOD,        ONLY : NTLOOP_PRIOR

#     include "CMN_SIZE"   ! Size parameters
#     include "comode.h"   ! NTLOOP, IGAS
#     include "CMN_VEL"    ! XYLAI

      ! Arguments
      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS

      ! Local Variables      
      INTEGER              :: I,    I0, IOS, J,  J0, L, N, JLOOP
      INTEGER              :: YYYY, MM, DD,  HH, SS, ZIP_HH
      INTEGER              :: IJLOOP
      CHARACTER(LEN=255)   :: FILENAME
      !>>>                
      ! Now include adjoint of F (dkh, 10/03/08) 
      INTEGER              :: NS
      !<<<


      ! Temporary storage arrays for checkpointed variables
      TYPE (XPLEX)               :: CHECK_RP_IN(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)               :: CHECK_FINAL(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)               :: CHECK_RP_OUT(IIPAR,JJPAR,LLPAR)
      ! Always recompute, so these don't need to be checkponted (dkh, 06/11/09) 
!      TYPE (XPLEX)               :: CHECK1(IIPAR,JJPAR,LLPAR)
!      TYPE (XPLEX)               :: CHECK2(IIPAR,JJPAR,LLPAR)
!      TYPE (XPLEX)               :: CHECK3(IIPAR,JJPAR,LLPAR)
!      TYPE (XPLEX)               :: CHECK4(IIPAR,JJPAR,LLPAR)
!      TYPE (XPLEX)               :: CHECK5(IIPAR,JJPAR,LLPAR)
!      TYPE (XPLEX)               :: CHECK6(IIPAR,JJPAR,LLPAR)
!      TYPE (XPLEX)               :: CHECK7(IIPAR,JJPAR,LLPAR)
!      TYPE (XPLEX)               :: CHECK8(IIPAR,JJPAR,LLPAR)
!      TYPE (XPLEX)               :: CHECK9(IIPAR,JJPAR,LLPAR)
      ! Now use NTLOOP because we want everything incase LVARTROP (dkh, 08/04/09) 
      !TYPE (XPLEX)               :: SMVGARRAY(NTT,IGAS)
      TYPE (XPLEX)               :: SMVGARRAY(NTLOOP,IGAS)

      ! For binary punch file, version 2.0
      TYPE (XPLEX)               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      INTEGER 		   :: MAX_nitr_max
      INTEGER 		   :: NSOFAR
      
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 


      !=================================================================
      ! MAKE_CHECKPT_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      OUTPUT_CHECKPT_FILE = 'gctm.chk.YYYYMMDD.hhmm'

      ! Clear some arrays 
      CHECK_RP_IN(:,:,:)  = 0e0
      CHECK_FINAL(:,:,:)  = 0e0
      CHECK_RP_OUT(:,:,:) = 0e0
      ! Always recompute, so these don't need to be checkponted (dkh, 06/11/09) 
!      CHECK1(:,:,:) 	  = 0e0
!      CHECK2(:,:,:) 	  = 0e0
!      CHECK3(:,:,:) 	  = 0e0
!      CHECK4(:,:,:) 	  = 0e0
!      CHECK5(:,:,:) 	  = 0e0
!      CHECK6(:,:,:) 	  = 0e0
!      CHECK7(:,:,:) 	  = 0e0
!      CHECK8(:,:,:) 	  = 0e0
!      CHECK9(:,:,:) 	  = 0e0
      SMVGARRAY(:,:) 	  = 0e0


      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM Checkpoint File: ' // 
     &           'Instantaneous Tracer Concentrations (v/v)'
      CATEGORY = 'IJ-CHK-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the checkpoint file for output -- binary punch format
      !=================================================================

      ! Copy the output checkpoint file name into a local variable
      FILENAME = TRIM( OUTPUT_CHECKPT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to filename
      FILENAME = TRIM( ADJTMP_DIR ) // TRIM( FILENAME )


      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_CHECKPT_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      !=================================================================
      ! Write each checkpointed quantity to the checkpoint file
      !=================================================================
    
      IF ( LSULF .and. LAERO_THERM ) THEN 
         ! First write the input to RPMARES
         DO N = 1, NRPIN
   
            ! Set UNIT
            IF ( N == 6 ) THEN
         
               ! RH
               UNIT = '%'  

            ELSEIF ( N == 7 ) THEN

               ! Temp
               UNIT = 'K' 
   
            ELSE
 
               ! Some concentration 
               UNIT = 'ug/m3'

            ENDIF           

 
            ! Temporarily store data in CHECK_RP_IN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               CHECK_RP_IN(I,J,L) = RP_IN(I,J,L,N)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,     JJPAR,     LLPAR,     I0+1,
     &                  J0+1,      1,         CHECK_RP_IN )
         ENDDO 

      ENDIF ! LSULF 

      ! Support for CH4 (kjw, dkh, 02/12/12, adj32_023) 
      IF ( .not. ITS_A_CH4_SIM() ) THEN

         ! Write the final concetration values as saved at the end of geos_mod.f
         UNIT = 'kg/box'
         ! Change to N_TRACERS (dkh, 06/11/09) 
         !DO N = 1, NOBS
         DO N = 1, N_TRACERS

            ! Temporarily store data in CHECK_FINAL
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               CHECK_FINAL(I,J,L) = CHK_STT(I,J,L,N)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N + NRPIN,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,     JJPAR,     LLPAR,     I0+1,
     &               J0+1,      1,         CHECK_FINAL )
         ENDDO

      ! It is a CH4 simulation 
      ELSE 
         ! Write the final concetration values as saved at the end of geos_mod.f
         UNIT = 'kg/box'
         DO N = 1, N_TRACERS
         
            ! Temporarily store data in CHECK_FINAL
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLPAR     
            DO J = 1, JJPAR     
            DO I = 1, IIPAR
               CHECK_FINAL(I,J,L) = STT(I,J,L,N)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            print*,'READ_CHECKPT: stt(14,14,14,1) =',
     &                    stt(14,14,14,1)
            print*,'READ_CHECKPT: check_final(14,14,14,1) =',
     &                    CHECK_FINAL(14,14,14)

            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  1,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,     JJPAR,     LLPAR,     I0+1,
     &                  J0+1,      1,         CHECK_FINAL )
         ENDDO
      ENDIF ! ITS_A_CH4_SIM  
      
      ! Checkpt additional values for full chem simulation
      ! Replace NSRCX (dkh, 06/11/09) 
      !IF ( NSRCX == 3 .AND. LCHEM ) THEN 
      IF ( ITS_A_FULLCHEM_SIM() .AND. LCHEM ) THEN 
 
         ! Write the final species concetrations after full chemistry
         UNIT = 'molec/cm3/box'
           
         ! Transfer to temp array so that we only checkpt NTLOOP values,
         ! not ITLOOP. 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( JLOOP, N )
         DO N     = 1, IGAS
         !DO JLOOP = 1, NTT
         DO JLOOP = 1, NTLOOP
          
            SMVGARRAY(JLOOP,N) = CSPEC_PRIOR(JLOOP,N)

         ENDDO
         ENDDO
!$OMP END PARALLEL DO
 
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  1 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               NTLOOP,    IGAS,      1,         I0+1,
     &               J0+1,      1,         SMVGARRAY )


         ! Set NSOFAR 
         NSOFAR = NSOFAR + 1

         ! Transfer to temp array so that we only checkpt NTLOOP values,
         ! not ITLOOP. 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( JLOOP )
         DO JLOOP = 1, NTT

            SMVGARRAY(JLOOP,1) = ( PART_CASE(JLOOP) )

         ENDDO
!$OMP END PARALLEL DO

         ! Checkpoint PART_CASE
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  1 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               NTT,       1,         1,             I0+1,
     &               J0+1,      1,         SMVGARRAY )


         ! Set NSOFAR 
         NSOFAR = NSOFAR + 1

         ! Write the tracer concetrations before chemisty
         UNIT = 'kg/box'
         DO N = 1, N_TRACERS

            ! Temporarily store data in CHECK_FINAL

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               CHECK_FINAL(I,J,L) = CHK_STT_BEFCHEM(I,J,L,N)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N + NSOFAR,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,     JJPAR,     LLPAR,     I0+1,
     &                  J0+1,      1,         CHECK_FINAL )
         ENDDO

         ! Set NSOFAR 
         NSOFAR = NSOFAR + N_TRACERS

         ! LVARTROP support for adj (dkh, 01/26/11)
         ! Write CSPEC_FULL_PRIOR 
         IF ( LVARTROP ) THEN
            UNIT = 'molec/cm3'
            DO N = 1, IGAS


               CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N + NSOFAR,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  ILONG,     ILAT,      IPVERT,    I0+1,
     &                  J0+1, 1,
     &          (CSPEC_FULL_PRIOR(1:ILONG,1:ILAT,1:IPVERT,N)))


            ENDDO

            ! Set NSOFAR 
            NSOFAR = NSOFAR + IGAS


            ! Write the 3-D to 1-D mappings 
            UNIT = 'none'
            CATEGORY = 'isave'

            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  1 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               NTLOOP_PRIOR,    3,       1,         I0+1,
     &               J0+1,      1,
     &          xplx(ISAVE_PRIOR(1:NTLOOP_PRIOR,:)) )


            ! Set NSOFAR 
            NSOFAR = NSOFAR + 1

            ! reset CATEGORY
            CATEGORY = 'IJ-CHK-$'


         ENDIF

         ! Write last internal time step used by Rosenbrock Solver
         ! (dkh, 09/06/05)  
         UNIT = 's'

         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  1 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     LLTROP,     I0+1,
     &               J0+1,      1,         CHK_HSAVE )


         ! Set NSOFAR 
         NSOFAR = NSOFAR + 1
      ENDIF 

      IF ( LSULF .AND. LCHEM ) THEN 
         ! Write the concentrations of SO2 and H2O2 used by CHEM_SO2
         ! (dkh, 10/12/05)
         UNIT = 'v/v'

         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  1 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     LLTROP,     I0+1,
     &               J0+1,      1,         SO2_CHK )

         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  2 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     LLTROP,     I0+1,
     &               J0+1,      1,         H2O2_CHK )

         ! Set NSOFAR 
         NSOFAR = NSOFAR + 2

      ENDIF 


      ! SOILNOX
      IF ( LSOILNOX ) THEN 
         UNIT = 'molec/cm2/s'
     
         ! Temporarily store data in CHECK_FINAL
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            CHECK_FINAL(I,J,1) = SOILNOX_CHK(I,J)
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! write to file
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  1 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     1,          I0+1,
     &               J0+1,      1,         CHECK_FINAL(:,:,1) )

         ! Update NSOFAR
         NSOFAR = NSOFAR + 1

      ENDIF

      ! Only do this for fullchem (mak, dkh, 01/06/10) 
      IF ( LCHEM .and. ITS_A_FULLCHEM_SIM() ) THEN 
 
         ! Now checkpoint XYLAI as well, as it is difficult to recalc
         DO N = 1, NTYPE 
 
            ! This mapping is clunky, but copied directly from rdlai.f 
            IJLOOP = 0
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               IJLOOP = IJLOOP + 1
               SMVGARRAY(IJLOOP,1) = ( XYLAI(IJLOOP,N) )
            END DO
            END DO

            ! Checkpoint XYLAI
            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  1 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               NTT,       1,         1,             I0+1,
     &               J0+1,      1,         SMVGARRAY )

            ! Set NSOFAR 
            NSOFAR = NSOFAR + 1

         ENDDO 


      ENDIF 


      ! SLBASE 
      IF ( LLIGHTNOX ) THEN 
         UNIT = 'molec/6h/box'
     
         ! Temporarily store data in CHECK_FINAL
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO l = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            CHECK_FINAL(I,J,L) = SLBASE_CHK(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! write to file
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  1 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     LLPAR,     I0+1,
     &               J0+1,      1,         CHECK_FINAL(:,:,:) )

         ! Update NSOFAR
         NSOFAR = NSOFAR + 1

      ENDIF

      ! Remove this, it wasn't a noticable improvement (dkh, 06/11/09) 
!      IF ( LADJ_TRAN ) THEN
!         UNIT = 'v/v'
!
!         ! CHK_STT_TD 
!         DO N = 1, NTRACE
!            ! Temporarily store data in CHECK_FINAL
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!            DO l = 1, LLPAR
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!               CHECK_FINAL(I,J,L) = CHK_STT_TD(I,J,L,N)
!            ENDDO
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!
!            ! write to file
!            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &                  HALFPOLAR, CENTER180, CATEGORY,  1 + NSOFAR,
!     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &                  IIPAR,     JJPAR,     LLPAR,     I0+1,
!     &                  J0+1,      1,         CHECK_FINAL(:,:,:) )
!
!            ! Update NSOFAR
!            NSOFAR = NSOFAR + 1
!         ENDDO
!         ! CHK_STT_TC 
!         DO N = 1, NTRACE
!            ! Temporarily store data in CHECK_FINAL
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!            DO l = 1, LLPAR
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!               CHECK_FINAL(I,J,L) = CHK_STT_TC(I,J,L,N)
!            ENDDO
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!            
!            ! write to file
!            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &                  HALFPOLAR, CENTER180, CATEGORY,  1 + NSOFAR,
!     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &                  IIPAR,     JJPAR,     LLPAR,     I0+1,
!     &                  J0+1,      1,         CHECK_FINAL(:,:,:) )
!            
!            ! Update NSOFAR 
!            NSOFAR = NSOFAR + 1
!         ENDDO
!      ENDIF

      ! Remove this obsolete option (which was always T) (dkh, 06/11/09) 
!      ! Check for recomputation -- if so, go ahead and finish up. 
!      IF (L_RECOMP) GOTO 444
!
!      ! It's been awhile since I've tried L_RECOMP = .FALSE. Some things
!      ! need to be update (NSOFAR, anything else?) dkh, 07/22/05
!      CALL ERROR_STOP( 'L_RECOMP = F not supported', 
!     &                 'MAKE_CHECKPT_FILE' ) 
!
!
!      ! Write the output from RPMARES 
!      DO N = 1, NRPOUT
!
!         ! Set UNIT
!         IF ( N == 9 ) THEN 
!  
!            ! EXIT value  
!            UNIT = 'unitless'
!
!         ELSE
! 
!            ! Some concentration
!            UNIT = 'ug/m3' 
!
!         ENDIF
!
!         ! Temporarily store data in CHECK_RPOUT
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!         DO L = 1, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            CHECK_RP_OUT(I,J,L) = RP_OUT(I,J,L,N)
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &               HALFPOLAR, CENTER180, CATEGORY,  N + NSOFAR,
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &               IIPAR,     JJPAR,     LLPAR,     I0+1,
!     &               J0+1,      1,         CHECK_RP_OUT )
!
!      ENDDO
!
!      NSOFAR = NSOFAR + NRPOUT
!
!      ! Write the values of nitr_max
!      CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &            HALFPOLAR, CENTER180, CATEGORY,  N + NSOFAR,
!     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &            IIPAR,     JJPAR,     LLPAR,     I0+1,
!     &            J0+1,      1,         REAL( nitr_max ) )
!
!      ! Calculate max of nitr_max
!      MAX_nitr_max = MAXVAL( nitr_max(:,:,:) ) 
!
!      ! Check to see that nitr_max is in the right range
!      IF ( MAX_nitr_max  > NNNMAX ) 
!     &   CALL ERROR_STOP( 'nitr_max > NNNMAX', 'MAKE_CHECKPT_FILE' ) 
!      IF ( MAX_nitr_max  == 0 )
!     &   CALL ERROR_STOP( 'MAXVAL (nitr_max) = 0', 'MAKE_CHECKPT_FILE' )
!
!      ! Now write the intermediate values necessary for adjoint computation
!      DO N = 1, MAX_nitr_max  
!  
!         ! Update tracer number   
!         NSOFAR  =  NRPIN + NOBS + 1 + NRPOUT + 1 + 9 * (N-1)
!         !Temporarily store quantities in the TRACER array
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!         DO L = 1, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            IF ( N <= nitr_max(I,J,L) ) THEN 
!               CHECK1(I,J,L) = gamaan_fwd (I,J,L,N) 
!               CHECK2(I,J,L) = gamold_fwd (I,J,L,N) 
!               CHECK3(I,J,L) = wh2o_fwd (I,J,L,N) 
!               CHECK4(I,J,L) = ynh4_fwd (I,J,L,N) 
!               CHECK5(I,J,L) = eror_fwd (I,J,L,N) 
!               CHECK6(I,J,L) = REAL ( exit_fwd (I,J,L,N) )
!               CHECK7(I,J,L) = gamana_fwd (I,J,L,N) 
!               CHECK8(I,J,L) = gamas1_fwd (I,J,L,N) 
!               CHECK9(I,J,L) = gamas2_fwd (I,J,L,N) 
!            ENDIF
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
!     &               HALFPOLAR, CENTER180, CATEGORY,  1 + NSOFAR,
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,   
!     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
!     &               J0+1,      1,         CHECK1 )
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
!     &               HALFPOLAR, CENTER180, CATEGORY,  2 + NSOFAR, 
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,   
!     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
!     &               J0+1,      1,         CHECK2 )
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
!     &               HALFPOLAR, CENTER180, CATEGORY,  3 + NSOFAR, 
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,   
!     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
!     &               J0+1,      1,         CHECK3 )
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
!     &               HALFPOLAR, CENTER180, CATEGORY,  4 + NSOFAR, 
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,   
!     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
!     &               J0+1,      1,         CHECK4 )
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
!     &               HALFPOLAR, CENTER180, CATEGORY,  5 + NSOFAR, 
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,   
!     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
!     &               J0+1,      1,         CHECK5 )
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
!     &               HALFPOLAR, CENTER180, CATEGORY,  6 + NSOFAR, 
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,   
!     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
!     &               J0+1,      1,         CHECK6 )
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
!     &               HALFPOLAR, CENTER180, CATEGORY,  7 + NSOFAR, 
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,   
!     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
!     &               J0+1,      1,         CHECK7 )
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
!     &               HALFPOLAR, CENTER180, CATEGORY,  8 + NSOFAR, 
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,   
!     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
!     &               J0+1,      1,         CHECK8 )
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,    
!     &               HALFPOLAR, CENTER180, CATEGORY,  9 + NSOFAR, 
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,   
!     &               IIPAR,     JJPAR,     LLPAR,     I0+1,            
!     &               J0+1,      1,         CHECK9 )
!      ENDDO  
!
 444  CONTINUE 

      ! Close file
      CLOSE( IU_RST )

      
       ! Remove obsolete option (dkh, 06/11/09) 
!      ! Zip files 
!      IF ( L_ZIP_CHECKPT ) CALL BATCH_ZIP( YYYYMMDD, HHMMSS, 'chk', 1 ) 
 
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_CHECKPT_FILE: wrote file' )

      ! Return to calling program
      END SUBROUTINE MAKE_CHECKPT_FILE

!------------------------------------------------------------------------------

      SUBROUTINE READ_CHECKPT_FILE( YYYYMMDD, HHMMSS ) 
!
!******************************************************************************
!  Subroutine READ_CHECKPT_FILE initializes GEOS-CHEM tracer concentrations 
!  from a checkpoint file (binary punch file format) 
!  (dkh, 8/30/04)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!
!  Passed via CMN:
!  ============================================================================
!  (1 ) TAU    : TAU value (elapsed hours) at start of diagnostic interval
!
!  Notes
!  (1 ) Just like READ_RESTART_FILE except
!	- load the variables from TRACER directly back into the CHECKPT array
!       - file name now includes hhmmss
!       - reads files from ADJ_DIR (and can unzip them if L_ZIP_CHECKPT) 
!       - removes .chk. files after reading (if L_DEL_CHECKPT)
!       dkh, 9/30/04
!  (2 ) Add DATE(2) and reference GET_NHMDe and GET_NHMSe to enable BATCH_ZIP
!        (dkh, 11/22/04)
!  (3 ) Add support for L_RECOMP option to recompute (rather than checkpoint)
!        variables RP_OUT etc.  (dkh, 02/09/05)
!  (4 ) Now read in values to CHK_STT (dkh, 03/03/05)
!  (5 ) Add CHK_PSC. (dkh, 03/16/05)
!  (6 ) Added support for full chemistry.  Add references to NVAR, CSPEC, JLOP, 
!       NTLOOP_FORKPP.
!       Add variables PART_CASE, JLOOP.  Disable L_RECOMP = FALSE option for now. 
!       Add SMVGARRAY
!       (dkh, 07/22/05)
!  (7 ) Add SO2_CHK and H2O2_CHK.  (dkh, 10/12/05)
!  (8 ) Add WETD_CHK_H2O2s_CHEMT, WETD_CHK_H2O2s_DYNT, etc. (dkh, 10/23/05) 
!  (9 ) Add WETD_CHK_SO2CHEMT, WETD_CHK_SO2_DYNT. (dkh, 10/31/05)  
!  (10) Add CONV_CHK_H2O2s_CHEMT, CONV_CHK_SO2s_CHEMT, etc.  (dkh, 11/22/05)  
!  (11) Add SOILNOX. (dkh, 02/06/07) 
!  (12) Move dynamic checkpointing to READ_CHK_DYN_FILE. (dkh, 02/01/09) 
!  (13) Remove obsolete options (L_DEL_CHECKPT, L_ZIP_CHECKPT, L_RECOMP), 
!        check for aeroosl simulation (LSULF) and update names to v8 (dkh, 06/11/09) 
!  (14) Add XYLAI (dkh, 10/14/09) 
!  (15) BUG FIX: LVARTROP treated correctly (dkh, 01/26/11)
!  (16) BUF FIX: Fill CSPEC with SMAL2 to prevent underflow later (dkh, 02/18/11) 
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,         ONLY : OPEN_BPCH2_FOR_READ
      USE COMODE_MOD,        ONLY : CHK_CSPEC, JLOP
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR 
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,          ONLY : IU_RST, IOERROR
      USE GCKPP_ADJ_GLOBAL,  ONLY : NTT
      USE LOGICAL_MOD,       ONLY : LCHEM , LSULF
      USE LOGICAL_MOD,       ONLY : LSOILNOX, LLIGHTNOX
      USE LOGICAL_MOD,       ONLY : LPRT
      USE LOGICAL_ADJ_MOD,   ONLY : LAERO_THERM
      USE LOGICAL_ADJ_MOD,   ONLY : LDEL_CHKPT
      USE RESTART_MOD,       ONLY : CHECK_DIMENSIONS
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TRACER_MOD,        ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,        ONLY : ITS_AN_AEROSOL_SIM
      USE TRACER_MOD,        ONLY : N_TRACERS
      USE UNIX_CMDS_MOD,     ONLY : REMOVE_CMD

      ! LVARTROP support for adj (dkh, 01/26/11)
      USE COMODE_MOD,        ONLY : CSPEC_FULL
      USE LOGICAL_MOD,       ONLY : LVARTROP
      USE COMODE_MOD,        ONLY : ISAVE_PRIOR


#     include "CMN_SIZE"   ! Size parameters
#     include "comode.h"   ! ITLOOP, IGAS
#     include "CMN_VEL"    ! XYLAI 

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER             :: I, IOS, J, L, N, JLOOP, NN, NTL
      INTEGER             :: NCOUNT(NNPAR) 
      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
      COMPLEX*16              :: D_TRACER(IIPAR,JJPAR,LLPAR)
      ! Remove these since we always recompute instead  
      ! of checkpointing (dkh, 06/11/09) 
!      TYPE (XPLEX)		  :: CHECK1(IIPAR,JJPAR,LLPAR) 
!      TYPE (XPLEX)		  :: CHECK2(IIPAR,JJPAR,LLPAR) 
!      TYPE (XPLEX)		  :: CHECK3(IIPAR,JJPAR,LLPAR) 
!      TYPE (XPLEX)		  :: CHECK4(IIPAR,JJPAR,LLPAR) 
!      TYPE (XPLEX)		  :: CHECK5(IIPAR,JJPAR,LLPAR) 
!      TYPE (XPLEX)		  :: CHECK6(IIPAR,JJPAR,LLPAR) 
!      TYPE (XPLEX)		  :: CHECK7(IIPAR,JJPAR,LLPAR) 
!      TYPE (XPLEX)		  :: CHECK8(IIPAR,JJPAR,LLPAR) 
!      TYPE (XPLEX)		  :: CHECK9(IIPAR,JJPAR,LLPAR) 
      TYPE (XPLEX)		  :: SMVGARRAY(ITLOOP,IGAS)
      COMPLEX*16                :: D_SMVGARRAY(ITLOOP,IGAS)
      !>>>                
      ! Now include adjoint of F (dkh, 10/03/08) 
      INTEGER             :: NS
      !<<<

      TYPE (XPLEX)              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME
      CHARACTER(LEN=255)  :: UNZIP_FILE_CMD
      CHARACTER(LEN=255)  :: REMOVE_CHK_FILE_CMD
     

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL,    NV
      INTEGER             :: IJLOOP
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      TYPE (XPLEX)              :: LONRES,    LATRES
      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
      COMPLEX*16              :: D_LONRES,    D_LATRES
      COMPLEX*16              :: D_ZTAU0,     D_ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED

      !=================================================================
      ! READ_CHECKPT_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      INPUT_CHECKPT_FILE = 'gctm.chk.YYYYMMDD.hhmm'

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0
      D_TRACER(:,:,:) = 0e0
      ! Remove these since we always recompute instead  
      ! of checkpointing (dkh, 06/11/09) 
!      CHECK1(:,:,:) = 0e0
!      CHECK2(:,:,:) = 0e0
!      CHECK3(:,:,:) = 0e0
!      CHECK4(:,:,:) = 0e0
!      CHECK5(:,:,:) = 0e0
!      CHECK6(:,:,:) = 0e0
!      CHECK7(:,:,:) = 0e0
!      CHECK8(:,:,:) = 0e0
!      CHECK9(:,:,:) = 0e0
      SMVGARRAY(:,:) = 0e0
      D_SMVGARRAY(:,:)=0e0
      !=================================================================
      ! Open checkpoint file and read top-of-file header
      !=================================================================
      
      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_CHECKPT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to name
      FILENAME = TRIM( ADJTMP_DIR ) // TRIM( FILENAME )


      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'C H E C K P T   F I L E   I N P U T'

      ! Obsolete (dkh, 06/11/09) 
!      ! Unzip checkpt file
!      IF ( L_ZIP_CHECKPT ) THEN
!         UNZIP_FILE_CMD = TRIM( GUNZIP_CMD ) // ' ' //    
!     &                  TRIM( FILENAME ) // ZIP_SUFFIX    
!         CALL SYSTEM( TRIM( UNZIP_FILE_CMD ) )
!         WRITE( 6, 99 ) TRIM( UNZIP_FILE_CMD ) 
! 99   FORMAT( '     - READ_CHECKPT_FILE: Executing: ',a )
!      ENDIF


      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_CHECKPT_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
      
      !=================================================================
      ! Read checkpointed variables 
      !=================================================================
 
      ! First read the input to RPMARES
      ! Add check for full chem with aerosols (dkh, 06/11/09) 
      IF ( LSULF .and. LAERO_THERM ) THEN 
         DO N = 1, NRPIN
           READ( IU_RST, IOSTAT=IOS ) 
     &        MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
           LONRES%r = dble(D_LONRES)
           LATRES%r = dble(D_LATRES)
           LONRES%i = dimag(D_LONRES)
           LATRES%i = dimag(D_LATRES)
            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) EXIT

            ! IOS > 0 is a real I/O error -- print error message
            IF ( IOS > 0 ) 
     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:4' )

            READ( IU_RST, IOSTAT=IOS ) 
     &           CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP
            ZTAU0%r = dble(D_ZTAU0)
           ZTAU1%r = dble(D_ZTAU1)
           ZTAU0%i = dimag(D_ZTAU0)
           ZTAU1%i = dimag(D_ZTAU1)
            IF ( IOS /= 0 ) 
     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:5')

            READ( IU_RST, IOSTAT=IOS ) 
     &           ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
          TRACER%r = dble(D_TRACER)
          TRACER%i = dimag(D_TRACER)
            IF ( IOS /= 0 ) 
     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:6')

        
            !==============================================================
            ! Assign data from the TRACER array to the STT array.
            !==============================================================
  
            ! Only process checkpoint data (i.e. mixing ratio)
            IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN 

               ! Make sure array dimensions are of global size
               ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
               CALL CHECK_DIMENSIONS( NI, NJ, NL )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  RP_IN(I,J,L,N) = TRACER(I,J,L)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

            ENDIF
         ENDDO

      ENDIF 

      ! Read the values of CHK_STT
      ! Change to N_TRACES (dkh, 06/11/09) 
      !DO N = 1, NOBS
      DO N = 1, N_TRACERS
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
       LONRES%r = dble(D_LONRES)
           LATRES%r = dble(D_LATRES)
           LONRES%i = dimag(D_LONRES)
           LATRES%i = dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:7' )

         READ( IU_RST, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0%r = dble(D_ZTAU0)
           ZTAU1%r = dble(D_ZTAU1)
           ZTAU0%i = dimag(D_ZTAU0)
           ZTAU1%i = dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:8' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
        TRACER%r = dble(D_TRACER)
          TRACER%i = dimag(D_TRACER)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:9' )

         !==============================================================
         ! Assign data from the TRACER array to the STT array.
         !==============================================================

         ! Only process checkpoint data (i.e. mixing ratio)
         IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN

            ! Make sure array dimensions are of global size
            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
            CALL CHECK_DIMENSIONS( NI, NJ, NL )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               CHK_STT(I,J,L,N) = TRACER(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ENDIF
      ENDDO

      ! Replace NSRCX (dkh, 06/11/09) 
      !IF ( NSRCX == 3 .AND. LCHEM ) THEN 
      IF ( ITS_A_FULLCHEM_SIM() .AND. LCHEM ) THEN 

         ! Read the values of CHK_CSPEC
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME,D_LONRES, D_LATRES, HALFPOLAR, CENTER180
        LONRES%r = dble(D_LONRES)
           LATRES%r = dble(D_LATRES)
           LONRES%i = dimag(D_LONRES)
           LATRES%i = dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 555

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:13' )

         READ( IU_RST, IOSTAT=IOS )
     &         CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &         NTL,      NN,       NL,   IFIRST, JFIRST, LFIRST,
     &         NSKIP
          ZTAU0%r = dble(D_ZTAU0)
           ZTAU1%r = dble(D_ZTAU1)
           ZTAU0%i = dimag(D_ZTAU0)
           ZTAU1%i = dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:14' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( D_SMVGARRAY(JLOOP,N), JLOOP=1,NTL ), N=1,NN ) 
          SMVGARRAY%r = dble(D_SMVGARRAY)
          SMVGARRAY%i = dimag(D_SMVGARRAY)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:16' )

         !==============================================================
         ! Assign data from the SMVGARRAY array to CHK_CSPEC
         !==============================================================

         ! Only process checkpoint data 
         IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN

            ! Check to make sure data is NTLOOPxNVAR 
            !Can't do this because RURALBOX hasn't been called for this
            !time step yet, so we don't know NTT yet (dkh, 07/31/09
            !CALL CHECK_DIMENSIONS_2( NTL,           NN,   NL, 
            !                         NTT,           IGAS, 1        )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( JLOOP, N )
            DO N     = 1, NN
            DO JLOOP = 1, NTLOOP
 
              ! BUG FIX: fill with 1d-99 to prevent underflow later (dkh, 02/18/11) 
              !CHK_CSPEC(JLOOP,N) = SMVGARRAY(JLOOP,N)
              CHK_CSPEC(JLOOP,N) = MAX(SMVGARRAY(JLOOP,N),SMAL2)

            ENDDO
            ENDDO
!$OMP END PARALLEL DO

! dkh debug
!            print*, ' In reac_checkpt: chk_cspec(FD) = ', 
!     &              CHK_CSPEC(JLOP(IFD,JFD,LFD),:)
!            print*, ' In reac_checkpt: smvgarray(FD) = ', 
!     &              SMVGARRAY(JLOP(IFD,JFD,LFD),:)
         ELSE
            CALL ERROR_STOP(' Category is not correct ', 
     &                      ' reading CHK_CSPEC, checkpt_mod')

         ENDIF

         ! Read in partition case PART_CASE
         READ( IU_RST, IOSTAT=IOS ) 
     &      MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
             LONRES%r = dble(D_LONRES)
           LATRES%r = dble(D_LATRES)
           LONRES%i = dimag(D_LONRES)
           LATRES%i = dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 555
            
         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 )
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:17' )
               
         READ( IU_RST, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NTL,          NN,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
            ZTAU0%r = dble(D_ZTAU0)
           ZTAU1%r = dble(D_ZTAU1)
           ZTAU0%i = dimag(D_ZTAU0)
           ZTAU1%i = dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:18' )
            
         READ( IU_RST, IOSTAT=IOS )
     &        ( ( D_SMVGARRAY(JLOOP,N), JLOOP=1,NTL ), N=1,NN ) 
          SMVGARRAY%r = dble(D_SMVGARRAY)
          SMVGARRAY%i = dimag(D_SMVGARRAY)
         IF ( IOS /= 0 )
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:20' )

         ! Convert from SMVGARRAY (REAL) to PART_CASE (INT)

!         ! Check to make sure data is NTLOOPx1 
!         CALL CHECK_DIMENSIONS_2( NTL,           NN,   NL,
!     &                            NTT            1,    1        )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( JLOOP )
         DO JLOOP = 1, NTL
               
            PART_CASE(JLOOP) = NINT( SMVGARRAY(JLOOP,NN) )

         ENDDO
!$OMP END PARALLEL DO

         ! Read the values of CHK_STT_BEFCHEM
         DO N = 1, N_TRACERS
            READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
          LONRES%r = dble(D_LONRES)
           LATRES%r = dble(D_LATRES)
           LONRES%i = dimag(D_LONRES)
           LATRES%i = dimag(D_LATRES)
            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) EXIT

            ! IOS > 0 is a real I/O error -- print error message
            IF ( IOS > 0 ) 
     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:21' )

            READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP
           ZTAU0%r = dble(D_ZTAU0)
           ZTAU1%r = dble(D_ZTAU1)
           ZTAU0%i = dimag(D_ZTAU0)
           ZTAU1%i = dimag(D_ZTAU1)
            IF ( IOS /= 0 ) 
     &         CALL IOERROR(IOS,IU_RST,'read_checkpt_file:22' )

            READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
            TRACER%r = dble(D_TRACER)
          TRACER%i = dimag(D_TRACER)
            IF ( IOS /= 0 ) 
     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:23' )

            !==============================================================
            ! Assign data from the TRACER array to the STT array.
            !==============================================================

            ! Only process checkpoint data (i.e. mixing ratio)
            IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN

               ! Make sure array dimensions are of global size
               ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
               CALL CHECK_DIMENSIONS( NI, NJ, NL )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  CHK_STT_BEFCHEM(I,J,L,N) = TRACER(I,J,L)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

            ENDIF
         ENDDO

         ! LVARTROP support for adj (dkh, 01/26/11)

         IF ( LVARTROP ) THEN 
         ! Read the values of CSPEC_FULL
         DO N = 1, IGAS
            READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
                LONRES%r = dble(D_LONRES)
           LATRES%r = dble(D_LATRES)
           LONRES%i = dimag(D_LONRES)
           LATRES%i = dimag(D_LATRES)
            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) EXIT

            ! IOS > 0 is a real I/O error -- print error message
            IF ( IOS > 0 )
     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:210' )
         
            READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP
           ZTAU0%r = dble(D_ZTAU0)
           ZTAU1%r = dble(D_ZTAU1)
           ZTAU0%i = dimag(D_ZTAU0)
           ZTAU1%i = dimag(D_ZTAU1)
            IF ( IOS /= 0 )
     &         CALL IOERROR(IOS,IU_RST,'read_checkpt_file:220' )
         
            READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
           TRACER%r = dble(D_TRACER)
          TRACER%i = dimag(D_TRACER)
            IF ( IOS /= 0 )
     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:230' )

            !==============================================================
            ! Assign data from the TRACER array to the CSPEC_FULL array.
            !==============================================================

            ! Only process checkpoint data (i.e. mixing ratio)
            IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN

               ! Make sure array dimensions are of global size
               ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
               !CALL CHECK_DIMENSIONS( NI, NJ, NL )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
               DO L = 1, IPVERT
               DO J = 1, ILAT
               DO I = 1, ILONG
                  ! BUG FIX: fill with 1d-99 to prevent underflow later (dkh, 02/18/11) 
                  !CSPEC_FULL(I,J,L,N) = TRACER(I,J,L)
                  CSPEC_FULL(I,J,L,N) = MAX(TRACER(I,J,L),SMAL2)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

            ENDIF
         ENDDO

         ! Read the values of ISAVE_PRIOR 
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LONRES%r = dble(D_LONRES)
           LATRES%r = dble(D_LATRES)
           LONRES%i = dimag(D_LONRES)
           LATRES%i = dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 555

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 )
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:13' )

         READ( IU_RST, IOSTAT=IOS )
     &         CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &         NTL,      NN,       NL,   IFIRST, JFIRST, LFIRST,
     &         NSKIP
         ZTAU0%r = dble(D_ZTAU0)
           ZTAU1%r = dble(D_ZTAU1)
           ZTAU0%i = dimag(D_ZTAU0)
           ZTAU1%i = dimag(D_ZTAU1)
         IF ( IOS /= 0 )
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:14' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( D_SMVGARRAY(JLOOP,N), JLOOP=1,NTL ), N=1,NN )
         SMVGARRAY%r = dble(D_SMVGARRAY)
         SMVGARRAY%i = dimag(D_SMVGARRAY)
         IF ( IOS /= 0 )
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:16' )

         !==============================================================
         ! Assign data from the SMVGARRAY array to CHK_CSPEC
         !==============================================================

         ! Only process checkpoint data 
         IF ( CATEGORY(1:8) == 'isave' ) THEN

            ! Check to make sure data is NTLOOPxNVAR 
            !Can't do this because RURALBOX hasn't been called for this
            !time step yet, so we don't know NTT yet (dkh, 07/31/09
            !CALL CHECK_DIMENSIONS_2( NTL,           NN,   NL, 
            !                         NTT,           IGAS, 1        )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( JLOOP, N )
            DO N     = 1, 3
            DO JLOOP = 1, NTL

               ISAVE_PRIOR(JLOOP,N) = NINT(SMVGARRAY(JLOOP,N))

            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE
            CALL ERROR_STOP(' Category is not correct ',
     &                      ' reading CHK_CSPEC, checkpt_mod')

         ENDIF

         ENDIF ! LVARTROP


         ! Read the values of CHK_HSAVE
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LONRES%r = dble(D_LONRES)
           LATRES%r = dble(D_LATRES)
           LONRES%i = dimag(D_LONRES)
           LATRES%i = dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 555

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:24' )

         READ( IU_RST, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0%r = dble(D_ZTAU0)
           ZTAU1%r = dble(D_ZTAU1)
           ZTAU0%i = dimag(D_ZTAU0)
           ZTAU1%i = dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:25' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( D_CHK_HSAVE(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
        CHK_HSAVE%r = dble(D_CHK_HSAVE)
        CHK_HSAVE%i = dimag(D_CHK_HSAVE)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:26' )
     
      ENDIF 

      IF ( LSULF .and. LCHEM ) THEN 
         ! Read the values of SO2_CHK
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
          LONRES%r = dble(D_LONRES)
           LATRES%r = dble(D_LATRES)
           LONRES%i = dimag(D_LONRES)
           LATRES%i = dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 555

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:27' )

         READ( IU_RST, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
          ZTAU0%r = dble(D_ZTAU0)
           ZTAU1%r = dble(D_ZTAU1)
           ZTAU0%i = dimag(D_ZTAU0)
           ZTAU1%i = dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:28' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( D_SO2_CHK(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
         SO2_CHK%r = dble(D_SO2_CHK)
         SO2_CHK%i = dimag(D_SO2_CHK)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:29' )

         ! Read the values of H2O2_CHK
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LONRES%r = dble(D_LONRES)
           LATRES%r = dble(D_LATRES)
           LONRES%i = dimag(D_LONRES)
           LATRES%i = dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 555

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:30' )

         READ( IU_RST, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
          ZTAU0%r = dble(D_ZTAU0)
           ZTAU1%r = dble(D_ZTAU1)
           ZTAU0%i = dimag(D_ZTAU0)
           ZTAU1%i = dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:31' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( D_H2O2_CHK(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
         H2O2_CHK%r = dble(D_H2O2_CHK)
         H2O2_CHK%i = dimag(D_H2O2_CHK)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:32' )

      ENDIF   ! LCHEM

      ! SOILNOX
      IF ( LSOILNOX ) THEN 
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LONRES%r = dble(D_LONRES)
           LATRES%r = dble(D_LATRES)
           LONRES%i = dimag(D_LONRES)
           LATRES%i = dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 555

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:63' )

         READ( IU_RST, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0%r = dble(D_ZTAU0)
           ZTAU1%r = dble(D_ZTAU1)
           ZTAU0%i = dimag(D_ZTAU0)
           ZTAU1%i = dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:64' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), 
     &              L=1,NL )
          TRACER%r = dble(D_TRACER)
          TRACER%i = dimag(D_TRACER)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:59' )

         SOILNOX_CHK(:,:) = TRACER(:,:,1)

      ENDIF 

      ! Only do this for fullchem (mak, dkh, 01/06/10) 
      IF ( LCHEM .and. ITS_A_FULLCHEM_SIM() ) THEN
      
         ! Read in partition case XYLAI (dkh, 10/14/09) 
         DO NV = 1 , NTYPE 

            READ( IU_RST, IOSTAT=IOS ) 
     &         MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
             LONRES%r = dble(D_LONRES)
           LATRES%r = dble(D_LATRES)
           LONRES%i = dimag(D_LONRES)
           LATRES%i = dimag(D_LATRES)
            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) GOTO 555
            
            ! IOS > 0 is a real I/O error -- print error message
            IF ( IOS > 0 )
     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:17' )
               
            READ( IU_RST, IOSTAT=IOS )
     &         CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &         NTL,          NN,       NL,   IFIRST, JFIRST, LFIRST,
     &         NSKIP
             ZTAU0%r = dble(D_ZTAU0)
           ZTAU1%r = dble(D_ZTAU1)
           ZTAU0%i = dimag(D_ZTAU0)
           ZTAU1%i = dimag(D_ZTAU1)
            IF ( IOS /= 0 )
     &         CALL IOERROR(IOS,IU_RST,'read_checkpt_file:18' )
          
            READ( IU_RST, IOSTAT=IOS )
     &           ( ( D_SMVGARRAY(JLOOP,N), JLOOP=1,NTL ), N=1,NN ) 
           SMVGARRAY%r = dble(D_SMVGARRAY)
          SMVGARRAY%i = dimag(D_SMVGARRAY)
            IF ( IOS /= 0 )
     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:20' )

            ! This mapping is clunky, but copied directly from rdlai.f 
            IJLOOP = 0
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               IJLOOP = IJLOOP + 1
               XYLAI(IJLOOP,NV) = SMVGARRAY(IJLOOP,1)
            END DO
            END DO

         ENDDO 
 
      ENDIF

      ! SLBASE
      IF ( LLIGHTNOX ) THEN 
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
          LONRES%r = dble(D_LONRES)
           LATRES%r = dble(D_LATRES)
           LONRES%i = dimag(D_LONRES)
           LATRES%i = dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 555

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:63' )

         READ( IU_RST, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0%r = dble(D_ZTAU0)
           ZTAU1%r = dble(D_ZTAU1)
           ZTAU0%i = dimag(D_ZTAU0)
           ZTAU1%i = dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:64' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), 
     &              L=1,NL )
          TRACER%r = dble(D_TRACER)
          TRACER%i = dimag(D_TRACER)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:59' )

         SLBASE_CHK(:,:,:) = TRACER(:,:,:)

      ENDIF 

      ! Take this part out, it didn't help much (dkh, 06/11/09) 
!      IF ( LADJ_TRAN ) THEN
!         ! Read the values of CHK_STT_TD
!         DO N = 1, NTRACE
!            READ( IU_RST, IOSTAT=IOS )
!     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!            ! IOS < 0 is end-of-file, so exit
!            IF ( IOS < 0 ) EXIT
!
!            ! IOS > 0 is a real I/O error -- print error message
!            IF ( IOS > 0 )
!     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:121' )
!
!            READ( IU_RST, IOSTAT=IOS )
!     &           CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &           NSKIP
!
!            IF ( IOS /= 0 )
!     &         CALL IOERROR(IOS,IU_RST,'read_checkpt_file:122' )
!
!            READ( IU_RST, IOSTAT=IOS )
!     &           ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!            IF ( IOS /= 0 )
!     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:123' )
!
!            !==============================================================
!            ! Assign data from the TRACER array to the STT array.
!            !==============================================================
!            
!            ! Only process checkpoint data (i.e. mixing ratio)
!            IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN
!               
!               ! Make sure array dimensions are of global size
!               ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
!               CALL CHECK_DIMENSIONS( NI, NJ, NL )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L ) 
!               DO L = 1, LLPAR
!               DO J = 1, JJPAR
!               DO I = 1, IIPAR
!                  CHK_STT_TD(I,J,L,N) = TRACER(I,J,L)
!               ENDDO
!               ENDDO
!               ENDDO
!!$OMP END PARALLEL DO
!            
!            ENDIF
!         ENDDO
!
!         ! Read the values of CHK_STT_TC
!         DO N = 1, NTRACE 
!            READ( IU_RST, IOSTAT=IOS )
!     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!            
!            ! IOS < 0 is end-of-file, so exit
!            IF ( IOS < 0 ) EXIT
!            
!            ! IOS > 0 is a real I/O error -- print error message
!            IF ( IOS > 0 ) 
!     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:221' )
!            
!            READ( IU_RST, IOSTAT=IOS )
!     &           CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &           NSKIP
!            
!            IF ( IOS /= 0 ) 
!     &         CALL IOERROR(IOS,IU_RST,'read_checkpt_file:222' )
!            
!            READ( IU_RST, IOSTAT=IOS )
!     &           ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!            
!            IF ( IOS /= 0 )  
!     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:223' )
!            
!            !==============================================================
!            ! Assign data from the TRACER array to the STT array.
!            !==============================================================
!            
!            ! Only process checkpoint data (i.e. mixing ratio)
!            IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN
!               
!               ! Make sure array dimensions are of global size
!               ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
!               CALL CHECK_DIMENSIONS( NI, NJ, NL )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L ) 
!               DO L = 1, LLPAR
!               DO J = 1, JJPAR
!               DO I = 1, IIPAR
!                  CHK_STT_TC(I,J,L,N) = TRACER(I,J,L)
!               ENDDO
!               ENDDO
!               ENDDO
!!$OMP END PARALLEL DO
!            
!            ENDIF
!         ENDDO
!      
!      ENDIF  !  LADJ_TRAN



      ! Take this out as I always had L_RECOMP = T 
!      ! Check for recomputation -- if so, go ahead and finish up
!      IF ( L_RECOMP ) GOTO 555
!
!      ! Read output from RPMARES  
!      DO N = 1, NRPOUT
!         READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!         ! IOS < 0 is end-of-file, so exit
!         IF ( IOS < 0 ) EXIT
!
!         ! IOS > 0 is a real I/O error -- print error message
!         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:4' )
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:5')
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:6')
!
!         !==============================================================
!         ! Assign data from the TRACER array to the STT array.
!         !==============================================================
! 
!         ! Only process checkpoint data (i.e. mixing ratio)
!         IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN
!
!            ! Make sure array dimensions are of global size
!            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
!            CALL CHECK_DIMENSIONS( NI, NJ, NL )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!            DO L = 1, LLPAR
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!               RP_OUT(I,J,L,N) = TRACER(I,J,L)
!            ENDDO
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!
!         ENDIF
!      ENDDO
!
!      ! Read nitr_max 
!      READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!      ! IOS < 0 is end-of-file, so return
!      IF ( IOS < 0 ) RETURN
!
!      ! IOS > 0 is a real I/O error -- print error message
!      IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:4' )
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:5')
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:6')
!
!      !==============================================================
!      ! Assign data from the TRACER array to the STT array.
!      !==============================================================
! 
!      ! Only process checkpoint data (i.e. mixing ratio)
!      IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN
!
!         ! Make sure array dimensions are of global size
!         ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
!         CALL CHECK_DIMENSIONS( NI, NJ, NL )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!         DO L = 1, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            nitr_max(I,J,L) = INT( TRACER(I,J,L) )
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!
!      ENDIF
!
!      ! Read the variables checkpointed for the adjoint calculation
!      ! CHECK1
!      READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!      ! IOS < 0 is end-of-file, so return
!      IF ( IOS < 0 ) RETURN
!
!      ! IOS > 0 is a real I/O error -- print error message
!      IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:4' )
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:5')
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( CHECK1(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:6')
!
!
!      ! CHECK2
!      READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!      ! IOS < 0 is end-of-file, so return
!      IF ( IOS < 0 ) RETURN
!
!      ! IOS > 0 is a real I/O error -- print error message
!      IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:4' )
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:5')
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( CHECK2(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:6')
!
!
!      ! CHECK3	
!      READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!      ! IOS < 0 is end-of-file, so return
!      IF ( IOS < 0 ) RETURN
!
!      ! IOS > 0 is a real I/O error -- print error message
!      IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:4' )
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:5')
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( CHECK3(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:6')
!         
!
!      ! CHECK4
!      READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!      ! IOS < 0 is end-of-file, so return
!      IF ( IOS < 0 ) RETURN
!
!      ! IOS > 0 is a real I/O error -- print error message
!      IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:4' )
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:5')
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( CHECK4(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:6')
!
!
!      ! CHECK5
!      READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!      ! IOS < 0 is end-of-file, so return
!      IF ( IOS < 0 ) RETURN
!
!      ! IOS > 0 is a real I/O error -- print error message
!      IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:4' )
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:5')
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( CHECK5(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:6')
!
!
!      ! CHECK6
!      READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!      ! IOS < 0 is end-of-file, so return
!      IF ( IOS < 0 ) RETURN
!
!      ! IOS > 0 is a real I/O error -- print error message
!      IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:4' )
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:5')
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( CHECK6(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:6')
!
!
!      ! CHECK7
!      READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!      ! IOS < 0 is end-of-file, so return
!      IF ( IOS < 0 ) RETURN
!
!      ! IOS > 0 is a real I/O error -- print error message
!      IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:4' )
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:5')
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( CHECK7(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:6')
! 
!
!      ! CHECK7 
!      READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!      ! IOS < 0 is end-of-file, so return
!      IF ( IOS < 0 ) RETURN
!
!      ! IOS > 0 is a real I/O error -- print error message
!      IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:4' )
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:5')
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( CHECK7(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:6')
!
!
!      ! CHECK8
!      READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!      ! IOS < 0 is end-of-file, so return
!      IF ( IOS < 0 ) RETURN
!
!      ! IOS > 0 is a real I/O error -- print error message
!      IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:4' )
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:5')
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( CHECK8(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:6')
!
!
!      ! CHECK9
!      READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!      ! IOS < 0 is end-of-file, so return
!      IF ( IOS < 0 ) RETURN
!
!      ! IOS > 0 is a real I/O error -- print error message
!      IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:4' )
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:5')
!
!      READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( CHECK9(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:6')
!
!
!      ! Write check arrays to the appropriate adjoint variables
!      DO N = 1, MAXVAL ( nitr_max(:,:,:) )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!         DO L = 1, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            gamaan_fwd (I,J,L,N) 	= CHECK1(I,J,L)
!            gamold_fwd (I,J,L,N) 	= CHECK2(I,J,L)
!            wh2o_fwd (I,J,L,N)		= CHECK3(I,J,L) 
!            ynh4_fwd (I,J,L,N)		= CHECK4(I,J,L) 
!            eror_fwd (I,J,L,N)		= CHECK5(I,J,L) 
!            exit_fwd (I,J,L,N)		= INT( CHECK6(I,J,L) )
!            gamana_fwd (I,J,L,N)	= CHECK7(I,J,L) 
!            gamas1_fwd (I,J,L,N)	= CHECK8(I,J,L) 
!            gamas2_fwd (I,J,L,N)	= CHECK9(I,J,L) 
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!
!      ENDDO
!
 555  CONTINUE

      ! Close file
      CLOSE( IU_RST )      

      ! Remove files if L_CHK_DEL = TRUE 
      IF ( LDEL_CHKPT ) THEN

        REMOVE_CHK_FILE_CMD  = TRIM ( REMOVE_CMD ) // ' ' //
     &                         TRIM ( FILENAME )

        CALL SYSTEM( TRIM( REMOVE_CHK_FILE_CMD ) )

        WRITE( 6, 102 ) TRIM( REMOVE_CHK_FILE_CMD )
 102    FORMAT( '     - READ_CHECKPT_FILE: Executing: ',a )

      ENDIF 

!      ! Zip the .chk. file if it hasn't been deleted and zipping 
!      ! is requested
!      IF ( L_ZIP_CHECKPT .AND. (.NOT. L_DEL_CHECKPT) ) THEN
!         CALL BATCH_ZIP( YYYYMMDD, HHMMSS, 'chk', -1 )
!      ENDIF


      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_CHECKPT_FILE: read file' )

      ! Return to calling program
      END SUBROUTINE READ_CHECKPT_FILE

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_OBS_FILE( YYYYMMDD, HHMMSS )
!
!******************************************************************************
!  Subroutine MAKE_OBS_FILE creates GEOS-CHEM observation files of tracer
!  mixing ratios (v/v) in binary punch file format.
!  (dkh, 9/01/04)
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a checkpoint file
!
!  Passed via CMN:
!  ============================================================================
!  (1 ) TAU    : TAU value (elapsed hours) at start of diagnostic interval
!
!  Passed via CMN_ADJ
!  ============================================================================
!  (1 ) CHECKPT : Array of quantities to be checkpointed    
!                               dim=(IIPAR,JJPAR,LLPAR,NCHECKPT)
!
!  NOTES:
!  (1 ) Just like MAKE_CHK_FILE except
!        - write to .obs. file
!        - only write output from rpmares,
!  (2 ) Switch to using OBS_STT rather than OBS
!  (3 ) Update to v8 format, remove obsolete options (dkh, 06/11/09)
!  (18) OBS_STT now in adj_arrays_mod.f instead of checkpt_mod.f (mak, 6/14/09)
!******************************************************************************
!    
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD,   ONLY : ADJTMP_DIR 
      USE ERROR_MOD,           ONLY : DEBUG_MSG
      USE FILE_MOD,            ONLY : IU_RST,      IOERROR
      USE GRID_MOD,            ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,         ONLY : LPRT 
      USE TIME_MOD,            ONLY : EXPAND_DATE, GET_TAU
      USE TRACER_MOD,          ONLY : N_TRACERS
      USE ADJ_ARRAYS_MOD,      ONLY : OBS_STT
      USE DIRECTORY_ADJ_MOD,   ONLY : ADJTMP_DIR

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS

      ! Local Variables     
      INTEGER              :: I,    I0, IOS, J,  J0, L, N
      INTEGER              :: YYYY, MM, DD,  HH, SS
      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      TYPE (XPLEX)               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT    
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE

      !=================================================================
      ! MAKE_OBS_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      OUTPUT_OBS_FILE = 'gctm.obs.YYYYMMDD.hhmm'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM OBS File: ' //
     &           'Observation Concentrations (kg/box)'
      UNIT     = 'kg/box'
      CATEGORY = 'IJ-OBS-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the observation file for output -- binary punch format
      !=================================================================

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_OBS_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to FILENAME
      FILENAME = TRIM( ADJTMP_DIR ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_OBS_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      !=================================================================
      ! Write each observed quantity to the observation file
      !=================================================================
 
      !DO N = 1, NOBS
      DO N = 1, N_TRACERS 

         !Temporarily store quantities in the TRACER array

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            TRACER(I,J,L) = OBS_STT(I,J,L,N)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,   
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,  
     &               IIPAR,     JJPAR,     LLPAR,     I0+1,
     &               J0+1,      1,         TRACER )

      ENDDO 

      ! Close file
      CLOSE( IU_RST )

      ! Obsolete (dkh, 06/11/09) 
!      ! Zip the obs file
!      IF ( L_ZIP_OBS ) THEN
!        
!         CALL BATCH_ZIP( YYYYMMDD, HHMMSS, 'obs', 1 )
!
!      ENDIF
 
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_OBS_FILE: wrote file' )

      ! Return to calling program
      END SUBROUTINE MAKE_OBS_FILE

!------------------------------------------------------------------------------

      SUBROUTINE READ_OBS_FILE( YYYYMMDD, HHMMSS )
!
!******************************************************************************
!  Subroutine READ_OBS_FILE reads the output of the reference run from an  
!  observation file (binary punch file format)
!  (dkh, 9/01/04)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!
!  Passed via CMN:
!  ============================================================================
!  (1 ) TAU    : TAU value (elapsed hours) at start of diagnostic interval
!
!  Notes
!  (1 ) Just like READ_CHECKPT_FILE except
!        - read NOBS variables into OBS array
!  (2 ) Switch to using OBS_STT rather than OBS (dkh 03/03/05)
!  (3 ) Update to v8 format, remove obsolete options (dkh, 06/11/09)
!  (4 ) OBS_STT now in adj_arrays_mod.f instead of checkpt_mod.f (mak, 6/14/09)
!******************************************************************************
!
      ! References to F90 modules
      USE RESTART_MOD,       ONLY : CHECK_DIMENSIONS
      USE BPCH2_MOD,         ONLY : OPEN_BPCH2_FOR_READ
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR 
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE FILE_MOD,          ONLY : IU_RST, IOERROR
      USE LOGICAL_MOD,       ONLY : LPRT
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TRACER_MOD,        ONLY : N_TRACERS
      USE ADJ_ARRAYS_MOD,    ONLY : OBS_STT

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER             :: I, IOS, J, L, N
      INTEGER             :: NCOUNT(NNPAR)
      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
      COMPLEX*16              :: D_TRACER(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME
      CHARACTER(LEN=255)  :: UNZIP_FILE_CMD
      CHARACTER(LEN=255)  :: ZIP_FILE_CMD

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      TYPE (XPLEX)              :: LONRES,    LATRES
      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
      COMPLEX*16              :: D_LONRES,    D_LATRES
      COMPLEX*16              :: D_ZTAU0,     D_ZTAU1

      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT    
      CHARACTER(LEN=40)   :: RESERVED

      !=================================================================
      ! READ_OBS_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      INPUT_OBS_FILE = 'gctm.obs.YYYYMMDD.hhmm'

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0
      D_TRACER(:,:,:)=0e0
      !=================================================================
      ! Open observation file and read top-of-file header
      !=================================================================

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_OBS_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to FILENAME
      !FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME ) 
      FILENAME = TRIM( ADJTMP_DIR ) // TRIM( FILENAME ) 

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'O B S   F I L E   I N P U T'

      ! Remove obsolete options  (dkh, 06/11/09) 
!      ! Unzip obs file
!      IF ( L_ZIP_OBS ) THEN
!         UNZIP_FILE_CMD = TRIM( GUNZIP_CMD ) // ' ' //
!     &                  TRIM( FILENAME ) // ZIP_SUFFIX
!         CALL SYSTEM( TRIM( UNZIP_FILE_CMD ) )
!         WRITE( 6, 99 ) TRIM( UNZIP_FILE_CMD )
! 99   FORMAT( '     - READ_OBS_FILE: Executing: ',a )
!      ENDIF

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_OBS_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )

      !=================================================================
      ! Read concentrations -- store in the TRACER array
      !=================================================================
      !DO N = 1, NOBS
      DO N = 1, N_TRACERS
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
          LONRES%r = dble(D_LONRES)
          LATRES%r = dble(D_LATRES)
          LONRES%i = dimag(D_LONRES)
          LATRES%i = dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_obs_file:4' )

         READ( IU_RST, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0%r = dble(D_ZTAU0)
          ZTAU1%r = dble(D_ZTAU1)
          ZTAU0%i = dimag(D_ZTAU0)
          ZTAU1%i = dimag(D_ZTAU1)
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_obs_file:5')

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
          TRACER%r = dble(D_TRACER)
          TRACER%i = dimag(D_TRACER)
         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_obs_file:6')

         !==============================================================
         ! Assign data from the TRACER array to the STT array.
         !==============================================================
 
         ! Only process observation data (i.e. aerosol and precursors)
         IF ( CATEGORY(1:8) == 'IJ-OBS-$' ) THEN

            ! Make sure array dimensions are of global size
            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
            CALL CHECK_DIMENSIONS( NI, NJ, NL )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               OBS_STT(I,J,L,N) = TRACER(I,J,L)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ENDIF
      ENDDO

      ! Close file
      CLOSE( IU_RST )     

      ! Remove obsolete options (dkh, 06/11/09) 
!      ! Zip the obs file
!      IF ( L_ZIP_OBS ) THEN
!
!         CALL BATCH_ZIP( YYYYMMDD, HHMMSS, 'obs', -1 ) 
!
!      ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_OBS_FILE: read file' )

      ! Return to calling program
      END SUBROUTINE READ_OBS_FILE

!----------------------------------------------------------------------- 
! Remove obsolete subroutines (dkh, 06/11/09) 
!      SUBROUTINE BATCH_ZIP( YYYYMMDD, HHMMSS, FID, MODE ) 
!!
!!**********************************************************************
!! Subroutine BATCH_ZIP zips a days worth of *.obs.*, *.adj.*, 
!!  and *.chk.* files using multiple processors.  Only works for
!!  TS_CHEM = 60 (min) and simulations that begin at HHMMSS = 000000. 
!!  Simulation ending at times other than HHMMSS = 000000 are allowed.
!!  The argument MODE indicates whether the batch of files to be zipped 
!!  begins (-1) or ends (+1) at HHMMSS, and adjustments are then made so
!!  that DATE(2) always indicates the time stamp of the latest file to be 
!!  zipped.   (dkh, 11/22/04)
!!
!! NOTES
!! 
!!**********************************************************************
!      ! Reference to f90 modules
!      USE TIME_MOD,   ONLY : GET_NYMDe,   GET_NHMSe,      EXPAND_DATE,
!     &                       GET_TS_CHEM, GET_TIME_AHEAD, GET_NHMSb
!      USE ERROR_MOD,  ONLY : ERROR_STOP
!
!#     include "CMN_SIZE"   !  Size parameters
!#     include "CMN_ADJ"    !  OBS_FREQ, GZIP_CMD 
!
!
!      ! Arguments
!      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS
!      CHARACTER(LEN=3)     :: FID
!      INTEGER		   :: MODE  	! +1 for fwd, -1 for backwd
!
!      ! Local variables
!      INTEGER		   :: HH, ZIP_INTERVAL, HH_MAX, ZHHMMSS
!      INTEGER		   :: NHMSe, NYMDe 
!      INTEGER              :: DATE(2) 
!      INTEGER 		   :: OBS_FREQ_HH = OBS_FREQ / 6 * 1000 
!      CHARACTER(LEN=255)   :: ZIP_FILE_CMD
!      CHARACTER(LEN=255)   :: ZIP_FILENAME
!      CHARACTER(LEN=255)   :: TO_ZIP_FILENAME
!
!      !------------------------------------------------------------
!      ! BATCH_ZIP begins here!
!      !------------------------------------------------------------ 
!
!      ! Check to make sure that TS_CHEM is actually 60 min 
!      ! and that the simulation began at the beginning of a day 
!      IF ( GET_TS_CHEM() /= 60 .OR. GET_NHMSb() /= 000000 ) THEN
!         WRITE(6,*) ' -- Timeing inappropriate for batch zip'
!         RETURN
!      ENDIF 
!
!      ! Get HHMMSS at end of run 
!      NHMSe = GET_NHMSe() 
!      NYMDe = GET_NYMDe()
!
!      ! Adjust the arguments YYYYMMDD and HHMMSS if we are operating
!      ! in reverse mode (i.e. zipping after reading)
!      IF ( MODE == -1 ) THEN
! 
!         ! Get YYYYMMDD and HHMMSS for 23 hours ahead
!         DATE = GET_TIME_AHEAD( 60 * 23 )
!
!         ! Adjust for case when is the zeroeth hour of final day
!         IF ( YYYYMMDD == NYMDe .AND. HHMMSS == 000000 ) THEN
!            
!            ! Set DATE(2) so that the final day's files get zipped
!            DATE(2) = NHMSe - 10000
!
!         ENDIF
!
!      ELSE
!        
!         DATE(1) = YYYYMMDD
!         DATE(2) = HHMMSS
!
!      ENDIF
! 
!      ! Determine range of batch of files to zip 
!      IF ( NHMSe == 000000 ) THEN
!
!         ! Batches will always span a full day
!         HH_MAX = 230000
!
!      ELSE 
!
!         ! Batch range depends upon the day 
!         IF ( YYYYMMDD == NYMDe ) THEN
!
!            ! The batch for the last day is shorter  
!            HH_MAX = NHMSe - 10000
!
!         ELSE 
!
!            ! Not the last day yet, so batch still spans a full day
!            HH_MAX = 230000
!
!         ENDIF
!
!      ENDIF 
!
!      IF ( FID == 'obs' ) THEN 
!         IF ( YYYYMMDD  == NYMDe .AND.
!     &      ( DATE(2) + OBS_FREQ_HH) > ( NHMSe - 10000 ) ) THEN
!            HH_MAX = DATE(2) 
!         ELSEIF ( ( DATE(2) + OBS_FREQ_HH ) > (230000) ) THEN
!            HH_MAX = DATE(2)
!         ENDIF 
!      ENDIF 
!
!      ! Only zip the batch of files at the end of the day (or partial day).
!      IF ( DATE(2) /= HH_MAX ) RETURN
!
!      ! Determine the number of files in the batch
!      IF ( FID == 'chk' .OR. FID == 'adj' ) THEN
!
!         ! There is (at least) one file for every hour
!         ZIP_INTERVAL = 10000
!
!      ELSEIF ( FID == 'obs' ) THEN
!
!            ! Convert the obseration interval (min) into the right units
!            ZIP_INTERVAL = OBS_FREQ_HH
!
!      ELSE
!        
!         CALL ERROR_STOP('File type not defined!', 
!     &                   'BATCH_ZIP (checkpt_mod.f)' )
!
!      ENDIF 
!
!      ! Create generic file name
!#if   defined( GEOS_1 ) || defined( GEOS_STRAT )
!      TO_ZIP_FILENAME = 'gctm.' // FID // '.YYYYMMDD.hhmmss'
!#else
!      TO_ZIP_FILENAME = 'gctm.' // FID // '.YYYYMMDD.hhmmss'
!#endif
! 
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( ZIP_FILENAME, HH      )
!!$OMP+PRIVATE( ZIP_FILE_CMD, ZHHMMSS )
!      DO HH = 000000, HH_MAX , ZIP_INTERVAL 
!
!         ! Set the HHMMSS of the file to be zipped
!         ZHHMMSS = HH 
!
!         ! Reconstruct name of file to be zipped
!         ZIP_FILENAME = TRIM(TO_ZIP_FILENAME)
! 
!         ! Replace YYYY, MM, DD, HH tokens in ZIP_FILENAME w/actual values
!         CALL EXPAND_DATE( ZIP_FILENAME, DATE(1), ZHHMMSS )
!
!         ! Add ADJ_DIR prefix to filename 
!         ZIP_FILENAME = TRIM( ADJ_DIR ) // TRIM( ZIP_FILENAME )
!
!         ! Create zip command
!         ZIP_FILE_CMD = TRIM( GZIP_CMD ) // ' ' //
!     &                  TRIM( ZIP_FILENAME )
!         CALL SYSTEM( TRIM ( ZIP_FILE_CMD ) )
!         WRITE( 6, 101 ) TRIM( ZIP_FILE_CMD )
! 101     FORMAT( '     - BATCH_ZIP: Executing: ',a )
!
!      ENDDO
!!$OMP END PARALLEL DO
!
!      ! Only continue when dealing with *.adj.* files
!      IF (FID /= 'adj' ) RETURN
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( ZIP_FILENAME, HH      )
!!$OMP+PRIVATE( ZIP_FILE_CMD, ZHHMMSS )
!      DO HH = 003000, HH_MAX + 3000, ZIP_INTERVAL 
!
!         ! Set the HHMMSS of the file to be zipped
!         ZHHMMSS = HH 
!
!         ! Replace YYYY, MM, DD, HH tokens in ZIP_FILENAME w/actual values
!         CALL EXPAND_DATE( ZIP_FILENAME, DATE(1), ZHHMMSS )
!   
!         ! Add ADJ_DIR prefix to filename 
!         ZIP_FILENAME = TRIM( ADJ_DIR ) // TRIM( ZIP_FILENAME )
!
!         ! Create zip command
!         ZIP_FILE_CMD = TRIM( GZIP_CMD ) // ' ' //
!     &                  TRIM( ZIP_FILENAME )
!         CALL SYSTEM( TRIM ( ZIP_FILE_CMD ) )
!         WRITE( 6, 101 ) TRIM( ZIP_FILE_CMD )
!
!      ENDDO
!!$OMP END PARALLEL DO
!
!
!      END SUBROUTINE BATCH_ZIP
!
!!----------------------------------------------------------------------

      SUBROUTINE CHECK_DIMENSIONS_2( XPASS, YPASS, ZPASS, 
     &                               XTRUE, YTRUE, ZTRUE )
!     
!******************************************************************************
!  Subroutine CHECK_DIMENSIONS_2makes sure that the dimensions of the
!  data for block that was checkpointed are correct.  XPASS should equal XTRUE,
!  etc.
!  (dkh, 07/22/05)
!
!  NOTES:
!  
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: XPASS, YPASS, ZPASS
      INTEGER, INTENT(IN) :: XTRUE, YTRUE, ZTRUE

      !=================================================================
      ! CHECK_DIMENSIONS_2 begins here!
      !=================================================================
      
      ! Error check longitude dimension: NI must equal IIPAR
      IF ( XPASS /= XTRUE .OR.
     &     YPASS /= YTRUE .OR.
     &     ZPASS /= ZTRUE      ) THEN  
         print*, XPASS, XTRUE 
         print*, YPASS, YTRUE 
         print*, ZPASS, ZTRUE 
         WRITE( 6, '(a)' ) 'ERROR reading in checkpt file!'
         WRITE( 6, '(a)' ) 'Wrong number of grid cells  encountered!'
         WRITE( 6, '(a)' ) 'STOP in CHECK_DIMENSIONS_2 (checkpt_mod.f)'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF
   
      ! Return to calling program
      END SUBROUTINE CHECK_DIMENSIONS_2

!----------------------------------------------------------------------
!
!      SUBROUTINE MAKE_SAVE_FILE( YYYYMMDD, HHMMSS, N_CALC )
!!
!!******************************************************************************
!!  Subroutine MAKE_SAVE_FILE creates GEOS-CHEM checkpt files of tracer 
!!  concentrations [kg/box].
!!  For use in checking chemistry adjoints.  (dkh, 07/19/06)  
!!  
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) YYYYMMDD : Year-Month-Date 
!!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a checkpoint file       
!!  (3 ) N_CALC   : Current iteration 
!!
!!  Passed via CMN:
!!  ============================================================================
!!  (1 ) TAU    : TAU value (elapsed hours) at start of diagnostic interval
!!
!!  Passed via ???:
!!  ============================================================================
!!  (1 ) CHECKPT	: Array of quantities to be checkpointed     
!!                               dim=(IIPAR,JJPAR,LLPAR,NCHECKPT) 
!!
!!  NOTES:
!!   Just like MAKE_CHECKPT_FILE except:
!!  
!!******************************************************************************
!!     
!      ! References to F90 modules
!      USE BPCH2_MOD
!      USE ERROR_MOD,  ONLY : DEBUG_MSG, ERROR_STOP
!      USE FILE_MOD,   ONLY : IU_RST,      IOERROR
!      USE GRID_MOD,   ONLY : GET_XOFFSET, GET_YOFFSET
!      USE TIME_MOD,   ONLY : EXPAND_DATE, GET_TAU
!      USE COMODE_MOD, ONLY : CSPEC_PRIOR , JLOP
!!      USE GCKPP_PARAMETERS, ONLY : NVAR
!      USE GCKPP_ADJ_GLOBAL,     ONLY : NTLOOP_FORKPP_ADJ
!
!#     include "CMN_SIZE"   ! Size parameters
!#     include "CMN"        ! TAU , NSRCX
!#     include "CMN_ADJ"    ! NRPIN, NRPOUT, L_ZIP_CHECKPT, GZIP_CMD, ADJ_DIR, NOBS
!#     include "CMN_SETUP"  ! LWETD
!#     include "comode.h"   ! IGAS
!
!      ! Arguments
!      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS, N_CALC
!
!      ! Local Variables      
!      INTEGER              :: I,    I0, IOS, J,  J0, L, N, JLOOP
!      INTEGER              :: YYYY, MM, DD,  HH, SS, ZIP_HH
!      CHARACTER(LEN=255)   :: FILENAME
!
!      ! Temporary storage arrays for checkpointed variables
!      TYPE (XPLEX)               :: CHECK_FINAL(IIPAR,JJPAR,LLPAR)
!
!      ! For binary punch file, version 2.0
!      TYPE (XPLEX)               :: LONRES, LATRES
!      INTEGER, PARAMETER   :: HALFPOLAR = 1
!      INTEGER, PARAMETER   :: CENTER180 = 1
!
!      CHARACTER(LEN=20)    :: MODELNAME
!      CHARACTER(LEN=40)    :: CATEGORY
!      CHARACTER(LEN=40)    :: UNIT     
!      CHARACTER(LEN=40)    :: RESERVED = ''
!      CHARACTER(LEN=80)    :: TITLE 
!
!
!      !=================================================================
!      ! MAKE_SAVE_FILE begins here!
!      !=================================================================
!
!      ! Hardwire output file for now
!#if   defined( GEOS_1 ) || defined( GEOS_STRAT )
!      OUTPUT_CHECKPT_FILE = 'gctm.save.YYMMDD.hhmmss.NN'
!#else
!      OUTPUT_CHECKPT_FILE = 'gctm.save.YYYYMMDD.hhmmss.NN'
!#endif
!
!      ! Clear some arrays 
!      CHECK_FINAL(:,:,:)  = 0e0
!
!      ! Define variables for BINARY PUNCH FILE OUTPUT
!      TITLE    = 'GEOS-CHEM Checkpoint File: ' // 
!     &           'Instantaneous Tracer Concentrations (v/v)'
!      CATEGORY = 'IJ-CHK-$'
!      LONRES   = DISIZE
!      LATRES   = DJSIZE
!
!      ! Call GET_MODELNAME to return the proper model name for
!      ! the given met data being used (bmy, 6/22/00)
!      MODELNAME = GET_MODELNAME()
!
!      ! Get the nested-grid offsets
!      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
!      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
!
!      !=================================================================
!      ! Open the save file for output -- binary punch format
!      !=================================================================
!
!      ! Copy the output checkpoint file name into a local variable
!      FILENAME = TRIM( OUTPUT_CHECKPT_FILE )
!
!      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
!
!      ! Append the iteration number suffix to the file name
!      CALL EXPAND_NAME( FILENAME, N_CALC )
!
!      ! Add ADJ_DIR prefix to filename
!      FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME )
!
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - MAKE_SAVE_FILE: Writing ', a )
!
!      ! Open checkpoint file for output
!      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
!
!      !=================================================================
!      ! Write each checkpointed quantity to the checkpoint file
!      !=================================================================
!    
!      ! Write the final concetration values as saved at the end of geos_mod.f
!      UNIT = 'kg/box'
!      DO N = 1, NOBS
!
!         ! Temporarily store data in CHECK_FINAL
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!         DO L = 1, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            CHECK_FINAL(I,J,L) = CHK_STT(I,J,L,N)
!         ENDDO
!         ENDDO
!         ENDDO
!!$OMP END PARALLEL DO
!
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &               HALFPOLAR, CENTER180, CATEGORY,  N,
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &               IIPAR,     JJPAR,     LLPAR,     I0+1,
!     &               J0+1,      1,         CHECK_FINAL )
!      ENDDO
!
!      ! Close file
!      CLOSE( IU_RST )
!      
!      ! Zip files 
!      IF ( L_ZIP_CHECKPT ) CALL BATCH_ZIP( YYYYMMDD, HHMMSS, 'save', 1 ) 
! 
!      !### Debug
!      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_SAVE_FILE: wrote file' )
!
!      ! Return to calling program
!      END SUBROUTINE MAKE_SAVE_FILE
!!-----------------------------------------------------------------------
!
!      SUBROUTINE MAKE_SAVE_FILE_2( YYYYMMDD, HHMMSS, N_CALC )
!!
!!******************************************************************************
!!  Subroutine MAKE_SAVE_FILE_2 creates GEOS-CHEM checkpt files of tracer 
!!  concentrations [kg/box]. Like MAKE_SAVE_FILE, except calculate the finite
!!  difference sensitivities directly.  Save these, the adjoint sensitivities, 
!!  and the ratio adj / fd.  Save first and 2nd order finite difference 
!!  sensitivities.  Requires running to XSTOP = 3. 
!!  For use in checking process specific adjoints. (dkh, 01/23/07) 
!!  
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) YYYYMMDD : Year-Month-Date 
!!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a checkpoint file       
!!  (3 ) N_CALC   : Current iteration 
!!
!!  Passed via CMN:
!!  ============================================================================
!!  (1 ) TAU    : TAU value (elapsed hours) at start of diagnostic interval
!!
!!  Passed via ???:
!!  ============================================================================
!!  (1 ) CHECKPT	: Array of quantities to be checkpointed     
!!                               dim=(IIPAR,JJPAR,LLPAR,NCHECKPT) 
!!
!!  NOTES:
!!   Just like MAKE_CHECKPT_FILE except:
!!  (1 ) Now write out both sets of 1st order FD gradients, and define a new 
!!        category (FD-TEST) for viewing in gamap.    (dkh, 10/10/08) 
!!******************************************************************************
!!     
!      ! References to F90 modules
!      USE BPCH2_MOD
!      USE ERROR_MOD,  ONLY : DEBUG_MSG, ERROR_STOP
!      USE FILE_MOD,   ONLY : IU_RST,      IOERROR
!      USE GRID_MOD,   ONLY : GET_XOFFSET, GET_YOFFSET
!      USE TIME_MOD,   ONLY : EXPAND_DATE, GET_TAU
!      USE COMODE_MOD, ONLY : CSPEC_PRIOR , JLOP
!!      USE GCKPP_PARAMETERS, ONLY : NVAR
!      USE GCKPP_ADJ_GLOBAL,     ONLY : NTLOOP_FORKPP_ADJ
!
!#     include "CMN_SIZE"   ! Size parameters
!#     include "CMN"        ! TAU , NSRCX
!#     include "CMN_ADJ"    ! NRPIN, NRPOUT, L_ZIP_CHECKPT, GZIP_CMD, ADJ_DIR, NOBS
!#     include "CMN_SETUP"  ! LWETD
!#     include "comode.h"   ! IGAS
!
!      ! Arguments
!      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS, N_CALC
!
!      ! Local Variables      
!      INTEGER              :: I,    I0, IOS, J,  J0, L, N, JLOOP
!      INTEGER              :: YYYY, MM, DD,  HH, SS, ZIP_HH
!      CHARACTER(LEN=255)   :: FILENAME
!      !CHARACTER(LEN=255)   :: INPUT_CHECKPT_FILE
!      CHARACTER(LEN=255)   :: INPUT_GDT_FILE
!      INTEGER              :: N_OFF(IIPAR,JJPAR)
!      TYPE (XPLEX), PARAMETER    :: FILTER = 1d0
!
!      ! Temporary storage arrays for checkpointed variables
!      TYPE (XPLEX)               :: CHK_STT_1(IIPAR,JJPAR,LLPAR,NOBS)
!      TYPE (XPLEX)               :: CHK_STT_2(IIPAR,JJPAR,LLPAR,NOBS)
!      TYPE (XPLEX)               :: CHK_STT_3(IIPAR,JJPAR,LLPAR,NOBS)
!      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
!      TYPE (XPLEX)               :: EMS_3D(IIPAR,JJPAR,MMSCL)
!      TYPE (XPLEX)               :: ADJ(IIPAR,JJPAR,1)
!      TYPE (XPLEX)               :: TRACER_2D(IIPAR,JJPAR,1)
!
!      ! For binary punch file, version 2.0
!      TYPE (XPLEX)               :: LONRES, LATRES
!      INTEGER              :: HALFPOLAR
!      INTEGER              :: CENTER180
!
!      CHARACTER(LEN=20)    :: MODELNAME
!      CHARACTER(LEN=40)    :: CATEGORY
!      CHARACTER(LEN=40)    :: UNIT     
!      CHARACTER(LEN=40)    :: RESERVED = ''
!      CHARACTER(LEN=80)    :: TITLE 
!
!      INTEGER             :: NI,     NJ,     NL
!      INTEGER             :: IFIRST, JFIRST, LFIRST
!      INTEGER             :: NTRACER,   NSKIP
!      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
!
!
!      !=================================================================
!      ! MAKE_SAVE_FILE_2 begins here!
!      !=================================================================
!
!
!      ! Clear some arrays 
!      CHK_STT_1(:,:,:,:)  = 0e0
!      CHK_STT_2(:,:,:,:)  = 0e0
!      CHK_STT_3(:,:,:,:)  = 0e0
!      EMS_3D(:,:,:)       = 0d0
!      ADJ(:,:,:)          = 0e0
!      TRACER(:,:,:)       = 0e0
!      TRACER_2D(:,:,:)    = 0e0
!      N_OFF(:,:)          = 0d0
!
!      !========================================
!      ! Read *.save* file from unperturbed run
!      !========================================
!      INPUT_CHECKPT_FILE = 'gctm.save.YYYYMMDD.hhmmss.NN'
!
!      ! Copy input file name to a local variable
!      FILENAME = TRIM( INPUT_CHECKPT_FILE )
!
!      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
!
!      ! Append the iteration number suffix to the file name
!      CALL EXPAND_NAME( FILENAME, 1 )
!
!      ! Add ADJ_DIR prefix to name
!      FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME )
!
!      ! Echo some input to the screen
!      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
!      WRITE( 6, '(a,/)' ) 'C H E C K P T   F I L E   I N P U T'
!
!
!      WRITE( 6, 400 ) TRIM( FILENAME )
! 400  FORMAT( '     - READ_CHECKPT_FILE: Reading ', a )
!
!      ! Open the binary punch file for input
!      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
!      
!      !=================================================================
!      ! Read checkpointed variables 
!      !=================================================================
!      ! Read the values of CHK_STT
!      DO N = 1, NOBS
!         READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!         ! IOS < 0 is end-of-file, so exit
!         IF ( IOS < 0 ) EXIT
!
!         ! IOS > 0 is a real I/O error -- print error message
!         IF ( IOS > 0 ) 
!     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:7' )
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!         IF ( IOS /= 0 ) 
!     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:8' )
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!         IF ( IOS /= 0 ) 
!     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:9' )
!
!         !==============================================================
!         ! Assign data from the TRACER array to the STT array.
!         !==============================================================
!
!         ! Only process checkpoint data (i.e. mixing ratio)
!         IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN
!
!            ! Make sure array dimensions are of global size
!            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
!            !CALL CHECK_DIMENSIONS( NI, NJ, NL )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!            DO L = 1, LLPAR
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!               CHK_STT_1(I,J,L,N) = TRACER(I,J,L)
!            ENDDO
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!
!         ENDIF
!      ENDDO
!      ! Close file
!      CLOSE( IU_RST )
!      
!
!      !========================================
!      ! Read *.save* file from perturbed run
!      !========================================
!      INPUT_CHECKPT_FILE = 'gctm.save.YYYYMMDD.hhmmss.NN'
!
!      ! Copy input file name to a local variable
!      FILENAME = TRIM( INPUT_CHECKPT_FILE )
!
!      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
!
!      ! Append the iteration number suffix to the file name
!      CALL EXPAND_NAME( FILENAME, 2 )
!
!      ! Add ADJ_DIR prefix to name
!      FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME )
!
!      ! Echo some input to the screen
!      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
!      WRITE( 6, '(a,/)' ) 'C H E C K P T   F I L E   I N P U T'
!
!
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - READ_CHECKPT_FILE: Reading ', a )
!
!      ! Open the binary punch file for input
!      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
!      
!      !=================================================================
!      ! Read checkpointed variables 
!      !=================================================================
!      ! Read the values of CHK_STT
!      DO N = 1, NOBS
!         READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!         ! IOS < 0 is end-of-file, so exit
!         IF ( IOS < 0 ) EXIT
!
!         ! IOS > 0 is a real I/O error -- print error message
!         IF ( IOS > 0 ) 
!     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:7' )
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!         IF ( IOS /= 0 ) 
!     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:8' )
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!         IF ( IOS /= 0 ) 
!     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:9' )
!
!         !==============================================================
!         ! Assign data from the TRACER array to the STT array.
!         !==============================================================
!
!         ! Only process checkpoint data (i.e. mixing ratio)
!         IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN
!
!            ! Make sure array dimensions are of global size
!            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
!            !CALL CHECK_DIMENSIONS( NI, NJ, NL )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!            DO L = 1, LLPAR
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!               CHK_STT_2(I,J,L,N) = TRACER(I,J,L)
!            ENDDO
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!
!         ENDIF
!      ENDDO
!      ! Close file
!      CLOSE( IU_RST )
!
!      !========================================
!      ! Read *.save* file from 2nd perturbed run
!      !========================================
!      INPUT_CHECKPT_FILE = 'gctm.save.YYYYMMDD.hhmmss.NN'
!
!      ! Copy input file name to a local variable
!      FILENAME = TRIM( INPUT_CHECKPT_FILE )
!
!      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
!
!      ! Append the iteration number suffix to the file name
!      CALL EXPAND_NAME( FILENAME, 3 )
!
!      ! Add ADJ_DIR prefix to name
!      FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME )
!
!      ! Echo some input to the screen
!      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
!      WRITE( 6, '(a,/)' ) 'C H E C K P T   F I L E   I N P U T'
!
!
!      WRITE( 6, 888 ) TRIM( FILENAME )
! 888  FORMAT( '     - READ_CHECKPT_FILE: Reading ', a )
!
!      ! Open the binary punch file for input
!      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
!      
!      !=================================================================
!      ! Read checkpointed variables 
!      !=================================================================
!      ! Read the values of CHK_STT
!      DO N = 1, NOBS
!         READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!         ! IOS < 0 is end-of-file, so exit
!         IF ( IOS < 0 ) EXIT
!
!         ! IOS > 0 is a real I/O error -- print error message
!         IF ( IOS > 0 ) 
!     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:7' )
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!         IF ( IOS /= 0 ) 
!     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:8' )
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!         IF ( IOS /= 0 ) 
!     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:9' )
!
!         !==============================================================
!         ! Assign data from the TRACER array to the STT array.
!         !==============================================================
!
!         ! Only process checkpoint data (i.e. mixing ratio)
!         IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN
!
!            ! Make sure array dimensions are of global size
!            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
!            !CALL CHECK_DIMENSIONS( NI, NJ, NL )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L )
!            DO L = 1, LLPAR
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!               CHK_STT_3(I,J,L,N) = TRACER(I,J,L)
!            ENDDO
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!
!         ENDIF
!      ENDDO
!      ! Close file
!      CLOSE( IU_RST )
!
!      !========================================
!      ! Read *.gdt.* file from unperturbed run
!      !========================================
!      ! Hardwire output file for now
!      INPUT_GDT_FILE = 'gctm.gdt.01'
!
!      ! Initialize some variables
!      TRACER(:,:,:) = 0e0
!
!      !=================================================================
!      ! Open gradient file and read top-of-file header
!      !=================================================================
!
!      ! Copy input file name to a local variable
!      FILENAME = TRIM( INPUT_GDT_FILE )
!
!      ! Add OPT_DATA_DIR prefix to FILENAME
!      FILENAME = TRIM( OPT_DATA_DIR ) // TRIM( FILENAME )
!
!      ! Echo some input to the screen
!      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
!      WRITE( 6, '(a,/)' ) 'G D T   F I L E   I N P U T'
!      WRITE( 6, 101 ) TRIM( FILENAME )
! 101  FORMAT( 'READ_GDT_FILE: Reading ', a )
!
!      ! Open the binary punch file for input
!      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
!
!         !=================================================================
!         ! Read adjoints -- store in the TRACER array
!         !=================================================================
!         DO N = 1, NNEMS
!            READ( IU_RST, IOSTAT=IOS )
!     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!            ! IOS < 0 is end-of-file, so exit
!            IF ( IOS < 0 ) EXIT
!
!            ! IOS > 0 is a real I/O error -- print error message
!            IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:5' )
!
!            READ( IU_RST, IOSTAT=IOS )
!     &           CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &           NSKIP
!
!            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:5')
!
!            READ( IU_RST, IOSTAT=IOS )
!     &           ( ( ( EMS_3D(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:7')
!
!            !==============================================================
!            ! Assign data from the TRACER array to the ADJ_STT array.
!            !==============================================================
!
!            ! Only process observation data (i.e. aerosol and precursors)
!            IF ( CATEGORY(1:8) == 'IJ-GDE-$' .and. N == EMSFD 
!     &           .and. MMSCL == 1 ) THEN
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!               DO J = 1, JJPAR
!               DO I = 1, IIPAR
!                  ADJ(I,J,1) = EMS_3D(I,J,1)
!               ENDDO
!               ENDDO
!!$OMP END PARALLEL DO
!
!            ENDIF
!         ENDDO
!      ! Close file
!      CLOSE( IU_RST )
!
!
!      !========================================
!      ! Write 2nd order FD gradient
!      !========================================
!      L = LFD
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         TRACER_2D(I,J,1) =(CHK_STT_2(I,J,L,NFD) - CHK_STT_3(I,J,L,NFD))
!     &          /  (2 * FD_DIFF )
!
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      ! Hardwire output file for now
!      OUTPUT_CHECKPT_FILE = 'gctm.save2.YYYYMMDD.hhmmss'
!
!      ! Define variables for BINARY PUNCH FILE OUTPUT
!      TITLE    = 'GEOS-CHEM Checkpoint File: ' // 
!     &           'Instantaneous Tracer Concentrations (v/v)'
!      CATEGORY = 'FD-TEST'
!      LONRES   = DISIZE
!      LATRES   = DJSIZE
!
!      ! Call GET_MODELNAME to return the proper model name for
!      ! the given met data being used (bmy, 6/22/00)
!      MODELNAME = GET_MODELNAME()
!
!      ! Get the nested-grid offsets
!      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
!      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
!
!      !=================================================================
!      ! Open the save file for output -- binary punch format
!      !=================================================================
!
!      ! Copy the output checkpoint file name into a local variable
!      FILENAME = TRIM( OUTPUT_CHECKPT_FILE )
!
!      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
!
!      ! Add ADJ_DIR prefix to filename
!      FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME )
!
!      WRITE( 6, 102 ) TRIM( FILENAME )
! 102  FORMAT( '     - MAKE_SAVE_FILE: Writing ', a )
!
!      ! Open checkpoint file for output
!      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
!
!      UNIT = 'kg/box'
!
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &               HALFPOLAR, CENTER180, CATEGORY,  1,
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &               IIPAR,     JJPAR,     1,         I0+1,
!     &               J0+1,      1,         TRACER_2D )
!
!      !========================================
!      ! Write ADJ gradient 
!      !========================================
!      UNIT = 'kg/box'
!
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &               HALFPOLAR, CENTER180, CATEGORY,  2,
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &               IIPAR,     JJPAR,     1,         I0+1,
!     &               J0+1,      1,         REAL(ADJ) )
!
!      !========================================
!      ! Write ADJ / FD ratio
!      !========================================
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         IF ( ABS(TRACER_2D(I,J,1)) .gt. FILTER ) THEN 
!
!            TRACER_2D(I,J,1) = REAL(ADJ(I,J,1)) / TRACER_2D(I,J,1)
!          
!         ELSE 
!
!            TRACER_2D(I,J,1) = 1d0
!  
!         ENDIF 
!
!         ! Keep track of number of points that are off
!         IF ( ( TRACER_2D(I,J,1) > 1d0 + FD_DIFF ) .OR.
!     &        ( TRACER_2D(I,J,1) < 1D0 - FD_DIFF )      ) THEN 
!            N_OFF(I,J) = 1
!         ENDIF 
!
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      UNIT = 'none'
!
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &               HALFPOLAR, CENTER180, CATEGORY,  3,
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &               IIPAR,     JJPAR,     1,         I0+1,
!     &               J0+1,      1,         TRACER_2D )
!
!      ! print out statistics of the ADJ / FD ratio 
!      WRITE(6,*) '===================================================='
!      WRITE(6,*) ' Global validation test for values > ', FILTER
!      WRITE(6,*) ' MAX of global 2nd order ADJ / FD  = ',
!     &         MAXVAL(TRACER_2D(:,:,1))
!      WRITE(6,*) ' MIN of global 2nd order ADJ / FD  = ',
!     &         MINVAL(TRACER_2D(:,:,1))
!      WRITE(6,*) ' Number of places where ratio off by ',FD_DIFF,' = ',
!     &         SUM(N_OFF(:,:))
!      WRITE(6,*) '===================================================='
!
!      !========================================
!      ! Write 1st order FD gradient
!      !========================================
!      UNIT = 'kg/box'
!      L = LFD
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         TRACER_2D(I,J,1) =(CHK_STT_2(I,J,L,NFD) - CHK_STT_1(I,J,L,NFD))
!     &          / FD_DIFF
!
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &            HALFPOLAR, CENTER180, CATEGORY,  4,
!     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &            IIPAR,     JJPAR,     1,         I0+1,
!     &            J0+1,      1,         TRACER_2D )
!
!      !========================================
!      ! Write chekpt values 
!      !========================================
!      UNIT = 'kg/box'
!      L = LFD
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         TRACER_2D(I,J,1) = CHK_STT_1(I,J,L,NFD)
!
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &            HALFPOLAR, CENTER180, CATEGORY,  5,
!     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &            IIPAR,     JJPAR,     1,         I0+1,
!     &            J0+1,      1,         TRACER_2D )
!
!      !========================================
!      ! Write the other 1st order FD gradient. (dkh, 10/10/08) 
!      !========================================
!      UNIT = 'kg/box'
!      L = LFD
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         TRACER_2D(I,J,1) =(CHK_STT_3(I,J,L,NFD) - CHK_STT_1(I,J,L,NFD))
!     &          / ( - FD_DIFF )
!
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &            HALFPOLAR, CENTER180, CATEGORY,  6,
!     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &            IIPAR,     JJPAR,     1,         I0+1,
!     &            J0+1,      1,         TRACER_2D )
!
!      ! Close file
!      CLOSE( IU_RST )
!      
!      ! Zip files 
!      IF ( L_ZIP_CHECKPT ) CALL BATCH_ZIP( YYYYMMDD, HHMMSS, 'save', 1 ) 
! 
!      !### Debug
!      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_SAVE_FILE_2: wrote file' )
!
!
!      ! Return to calling program
!      END SUBROUTINE MAKE_SAVE_FILE_2
!!----------------------------------------------------------------------
!
!      SUBROUTINE MAKE_SAVE_FILE_3( YYYYMMDD, N_CALC )
!!
!!******************************************************************************
!!  Subroutine MAKE_SAVE_FILE_3 creates GEOS-CHEM checkpt files of radiative
!!  forcing [kg/box]. Like MAKE_SAVE_FILE_2, it saves these, the adjoint 
!!  sensitivities,  and the ratio adj / fd.  Requires running to XSTOP = 3. 
!!  For use in checking process specific adjoints. (dkh, 07/09/08) 
!!  
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) YYYYMMDD : Year-Month-Date 
!!  (2 ) N_CALC   : Current iteration 
!!
!!  Passed via CMN:
!!  ============================================================================
!!  (1 ) TAU    : TAU value (elapsed hours) at start of diagnostic interval
!!
!!  Passed via ???:
!!  ============================================================================
!!  (1 ) CHECKPT	: Array of quantities to be checkpointed     
!!                               dim=(IIPAR,JJPAR,LLPAR,NCHECKPT) 
!!
!!  NOTES:
!!   Just like MAKE_CHECKPT_FILE except:
!!  
!!******************************************************************************
!!     
!      ! References to F90 modules
!      USE BPCH2_MOD
!      USE ERROR_MOD,  ONLY : DEBUG_MSG, ERROR_STOP
!      USE FILE_MOD,   ONLY : IU_RST,      IOERROR
!      USE GRID_MOD,   ONLY : GET_XOFFSET, GET_YOFFSET
!      USE TIME_MOD,   ONLY : EXPAND_DATE, GET_TAU
!      USE COMODE_MOD, ONLY : CSPEC_PRIOR , JLOP
!!      USE GCKPP_PARAMETERS, ONLY : NVAR
!      USE GCKPP_ADJ_GLOBAL,     ONLY : NTLOOP_FORKPP_ADJ
!
!#     include "CMN_SIZE"   ! Size parameters
!#     include "CMN"        ! TAU , NSRCX
!#     include "CMN_ADJ"    ! NRPIN, NRPOUT, L_ZIP_CHECKPT, GZIP_CMD, ADJ_DIR, NOBS
!#     include "CMN_SETUP"  ! LWETD
!#     include "comode.h"   ! IGAS
!
!      ! Arguments
!      INTEGER, INTENT(IN)  :: YYYYMMDD, N_CALC
!
!      ! Local Variables      
!      INTEGER              :: I,    I0, IOS, J,  J0, L, N, W, JLOOP
!      INTEGER              :: YYYY, MM, DD,  HH, SS, ZIP_HH
!      INTEGER              :: HHMMSS_dum
!      CHARACTER(LEN=255)   :: FILENAME
!      CHARACTER(LEN=255)   :: INPUT_GDT_FILE
!      CHARACTER(LEN=255)   :: INPUT_AOD_FILE
!      INTEGER              :: N_OFF(IIPAR,JJPAR)
!      TYPE (XPLEX), PARAMETER    :: FILTER = 1d-10
!
!      ! Temporary storage arrays for checkpointed variables
!      TYPE (XPLEX)               :: CHK_RAD_1(IIPAR,JJPAR)
!      TYPE (XPLEX)               :: CHK_RAD_2(IIPAR,JJPAR)
!      TYPE (XPLEX)               :: CHK_RAD_3(IIPAR,JJPAR)
!      TYPE (XPLEX)               :: TRACER(IIPAR,JJPAR,LLPAR)
!      TYPE (XPLEX)               :: EMS_3D(IIPAR,JJPAR,MMSCL)
!      TYPE (XPLEX)               :: ADJ(IIPAR,JJPAR,1)
!      TYPE (XPLEX)               :: TRACER_2D(IIPAR,JJPAR,1)
!
!      ! For binary punch file, version 2.0
!      TYPE (XPLEX)               :: LONRES, LATRES
!      INTEGER              :: HALFPOLAR
!      INTEGER              :: CENTER180
!
!      CHARACTER(LEN=20)    :: MODELNAME
!      CHARACTER(LEN=40)    :: CATEGORY
!      CHARACTER(LEN=40)    :: UNIT     
!      CHARACTER(LEN=40)    :: RESERVED = ''
!      CHARACTER(LEN=80)    :: TITLE 
!
!      INTEGER             :: NI,     NJ,     NL
!      INTEGER             :: IFIRST, JFIRST, LFIRST
!      INTEGER             :: NTRACER,   NSKIP
!      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
!
!      INTEGER, PARAMETER  :: NWL_MAX = 5 
!
!      !=================================================================
!      ! MAKE_SAVE_FILE_3 begins here!
!      !=================================================================
!
!
!      ! Clear some arrays 
!      CHK_RAD_1(:,:)      = 0e0
!      CHK_RAD_2(:,:)      = 0e0
!      CHK_RAD_3(:,:)      = 0e0
!      EMS_3D(:,:,:)       = 0d0
!      ADJ(:,:,:)          = 0e0
!      TRACER(:,:,:)       = 0e0
!      TRACER_2D(:,:,:)    = 0e0
!      N_OFF(:,:)          = 0d0
!
!      !========================================
!      ! Read *.save* file from unperturbed run
!      !========================================
!      !INPUT_CHECKPT_FILE = 'gctm.save.YYYYMMDD.hhmmss.NN'
!      INPUT_AOD_FILE = 'aod.YYYYMMDD.NN'
!
!      ! Copy input file name to a local variable
!      FILENAME = TRIM( INPUT_AOD_FILE )
!
!      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS_dum )
!
!      ! Append the iteration number suffix to the file name
!      CALL EXPAND_NAME( FILENAME, 1 )
!
!      ! Add ADJ_DIR prefix to name
!      FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME )
!
!      ! Echo some input to the screen
!      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
!      WRITE( 6, '(a,/)' ) 'A O D    F I L E   I N P U T'
!
!
!      WRITE( 6, 401 ) TRIM( FILENAME )
! 401  FORMAT( '     - MAKE_SAVE_FILE_3: Reading ', a )
!
!      ! Open the binary punch file for input
!      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
!      
!      !=================================================================
!      ! Read checkpointed variables 
!      !=================================================================
!      ! Read the values of CHK_STT
!      DO W = 1, NWL_MAX * 3 + 1
!         READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!         ! IOS < 0 is end-of-file, so exit
!         IF ( IOS < 0 ) EXIT
!
!         ! IOS > 0 is a real I/O error -- print error message
!         IF ( IOS > 0 ) 
!     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:7' )
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!         IF ( IOS /= 0 ) 
!     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:8' )
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!         IF ( IOS /= 0 ) 
!     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:9' )
!
!         !==============================================================
!         ! Assign data from the TRACER array to the STT array.
!         !==============================================================
!
!         ! Only process checkpoint data (i.e. mixing ratio)
!         IF ( CATEGORY(1:8) == 'IJ-AOD-$' .and.
!     &                    W == NWL_MAX * 3 + 1  ) THEN 
!
!            ! Make sure array dimensions are of global size
!            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
!            !CALL CHECK_DIMENSIONS( NI, NJ, NL )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!               CHK_RAD_1(I,J) = TRACER(I,J,1)
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!
!         ENDIF
!      ENDDO
!      ! Close file
!      CLOSE( IU_RST )
!      
!
!      !========================================
!      ! Read *.save* file from perturbed run
!      !========================================
!      !INPUT_CHECKPT_FILE = 'gctm.save.YYYYMMDD.hhmmss.NN'
!      INPUT_AOD_FILE = 'aod.YYYYMMDD.NN'
!
!      ! Copy input file name to a local variable
!      FILENAME = TRIM( INPUT_AOD_FILE )
!
!      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS_dum )
!
!      ! Append the iteration number suffix to the file name
!      CALL EXPAND_NAME( FILENAME, 2 )
!
!      ! Add ADJ_DIR prefix to name
!      FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME )
!
!      ! Echo some input to the screen
!      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
!      WRITE( 6, '(a,/)' ) 'A O D    F I L E   I N P U T'
!
!
!      WRITE( 6, 101 ) TRIM( FILENAME )
! 101  FORMAT( '     - MAKE_SAVE_FILE_3: Reading ', a )
!
!      ! Open the binary punch file for input
!      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
!      
!      !=================================================================
!      ! Read checkpointed variables 
!      !=================================================================
!      ! Read the values of CHK_STT
!      DO W = 1, NWL_MAX*3 + 1 
!         READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!         ! IOS < 0 is end-of-file, so exit
!         IF ( IOS < 0 ) EXIT
!
!         ! IOS > 0 is a real I/O error -- print error message
!         IF ( IOS > 0 ) 
!     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:7' )
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!         IF ( IOS /= 0 ) 
!     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:8' )
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!         IF ( IOS /= 0 ) 
!     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:9' )
!
!         !==============================================================
!         ! Assign data from the TRACER array to the STT array.
!         !==============================================================
!
!         ! Only process checkpoint data (i.e. mixing ratio)
!         IF ( CATEGORY(1:8) == 'IJ-AOD-$' .and. 
!     &                    W == NWL_MAX * 3 + 1  ) THEN
!
!            ! Make sure array dimensions are of global size
!            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
!            !CALL CHECK_DIMENSIONS( NI, NJ, NL )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!               CHK_RAD_2(I,J) = TRACER(I,J,1)
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!
!         ENDIF
!      ENDDO
!      ! Close file
!      CLOSE( IU_RST )
!
!      !========================================
!      ! Read *.save* file from 2nd perturbed run
!      !========================================
!      !INPUT_CHECKPT_FILE = 'gctm.save.YYYYMMDD.hhmmss.NN'
!      INPUT_AOD_FILE = 'aod.YYYYMMDD.NN'
!
!      ! Copy input file name to a local variable
!      FILENAME = TRIM( INPUT_AOD_FILE )
!
!      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS_dum )
!
!      ! Append the iteration number suffix to the file name
!      CALL EXPAND_NAME( FILENAME, 3 )
!
!      ! Add ADJ_DIR prefix to name
!      FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME )
!
!      ! Echo some input to the screen
!      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
!      WRITE( 6, '(a,/)' ) 'A O D   F I L E   I N P U T'
!
!
!      WRITE( 6, 889 ) TRIM( FILENAME )
! 889  FORMAT( '     - MAKE_SAVE_FILE_3: Reading ', a )
!
!      ! Open the binary punch file for input
!      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
!      
!      !=================================================================
!      ! Read checkpointed variables 
!      !=================================================================
!      ! Read the values of CHK_STT
!      DO W = 1, NWL_MAX * 3 + 1 
!         READ( IU_RST, IOSTAT=IOS )
!     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!         ! IOS < 0 is end-of-file, so exit
!         IF ( IOS < 0 ) EXIT
!
!         ! IOS > 0 is a real I/O error -- print error message
!         IF ( IOS > 0 ) 
!     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:7' )
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &        NSKIP
!
!         IF ( IOS /= 0 ) 
!     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:8' )
!
!         READ( IU_RST, IOSTAT=IOS )
!     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!         IF ( IOS /= 0 ) 
!     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:9' )
!
!         !==============================================================
!         ! Assign data from the TRACER array to the STT array.
!         !==============================================================
!
!         ! Only process checkpoint data (i.e. mixing ratio)
!         IF ( CATEGORY(1:8) == 'IJ-AOD-$' .and. 
!     &                    W == NWL_MAX * 3 + 1  ) THEN
!
!            ! Make sure array dimensions are of global size
!            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
!            !CALL CHECK_DIMENSIONS( NI, NJ, NL )
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!               CHK_RAD_3(I,J) = TRACER(I,J,1)
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!
!         ENDIF
!      ENDDO
!      ! Close file
!      CLOSE( IU_RST )
!
!      !========================================
!      ! Read *.gdt.* file from unperturbed run
!      !========================================
!      ! Hardwire output file for now
!      INPUT_GDT_FILE = 'gctm.gdt.01'
!
!      ! Initialize some variables
!      TRACER(:,:,:) = 0e0
!
!      !=================================================================
!      ! Open gradient file and read top-of-file header
!      !=================================================================
!
!      ! Copy input file name to a local variable
!      FILENAME = TRIM( INPUT_GDT_FILE )
!
!      ! Add OPT_DATA_DIR prefix to FILENAME
!      FILENAME = TRIM( OPT_DATA_DIR ) // TRIM( FILENAME )
!
!      ! Echo some input to the screen
!      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
!      WRITE( 6, '(a,/)' ) 'G D T   F I L E   I N P U T'
!      WRITE( 6, 109 ) TRIM( FILENAME )
! 109  FORMAT( 'READ_GDT_FILE: Reading ', a )
!
!      ! Open the binary punch file for input
!      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
!
!         !=================================================================
!         ! Read adjoints -- store in the TRACER array
!         !=================================================================
!         DO N = 1, NNEMS
!            READ( IU_RST, IOSTAT=IOS )
!     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
!
!            ! IOS < 0 is end-of-file, so exit
!            IF ( IOS < 0 ) EXIT
!
!            ! IOS > 0 is a real I/O error -- print error message
!            IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:5' )
!
!            READ( IU_RST, IOSTAT=IOS )
!     &           CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
!     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
!     &           NSKIP
!
!            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:5')
!
!            READ( IU_RST, IOSTAT=IOS )
!     &           ( ( ( EMS_3D(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
!
!            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:7')
!
!            !==============================================================
!            ! Assign data from the TRACER array to the ADJ_STT array.
!            !==============================================================
!
!            ! Only process observation data (i.e. aerosol and precursors)
!            IF ( CATEGORY(1:8) == 'IJ-GDE-$' .and. N == EMSFD 
!     &           .and. MMSCL == 1 ) THEN
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!               DO J = 1, JJPAR
!               DO I = 1, IIPAR
!                  ADJ(I,J,1) = EMS_3D(I,J,1)
!               ENDDO
!               ENDDO
!!$OMP END PARALLEL DO
!
!            ENDIF
!         ENDDO
!      ! Close file
!      CLOSE( IU_RST )
!
!
!      !========================================
!      ! Write 2nd order FD gradient
!      !========================================
!      !L = LFD
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         TRACER_2D(I,J,1) =(CHK_RAD_2(I,J) - CHK_RAD_3(I,J))
!     &          /  (2 * FD_DIFF )
!
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      ! Hardwire output file for now
!      OUTPUT_CHECKPT_FILE = 'gctm.save2.YYYYMMDD'
!
!      ! Define variables for BINARY PUNCH FILE OUTPUT
!      TITLE    = 'GEOS-CHEM Checkpoint File: ' // 
!     &           'Instantaneous Tracer Concentrations (v/v)'
!      CATEGORY = 'IJ-CHK-$'
!      LONRES   = DISIZE
!      LATRES   = DJSIZE
!
!      ! Call GET_MODELNAME to return the proper model name for
!      ! the given met data being used (bmy, 6/22/00)
!      MODELNAME = GET_MODELNAME()
!
!      ! Get the nested-grid offsets
!      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
!      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
!
!      !=================================================================
!      ! Open the save file for output -- binary punch format
!      !=================================================================
!
!      ! Copy the output checkpoint file name into a local variable
!      FILENAME = TRIM( OUTPUT_CHECKPT_FILE )
!
!      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
!      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS_dum )
!
!      ! Add ADJ_DIR prefix to filename
!      FILENAME = TRIM( ADJ_DIR ) // TRIM( FILENAME )
!
!      WRITE( 6, 102 ) TRIM( FILENAME )
! 102  FORMAT( '     - MAKE_SAVE_FILE: Writing ', a )
!
!      ! Open checkpoint file for output
!      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
!
!      UNIT = 'kg/box'
!
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &               HALFPOLAR, CENTER180, CATEGORY,  1,
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &               IIPAR,     JJPAR,     1,         I0+1,
!     &               J0+1,      1,         TRACER_2D )
!
!      !========================================
!      ! Write ADJ gradient 
!      !========================================
!      UNIT = 'kg/box'
!
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &               HALFPOLAR, CENTER180, CATEGORY,  2,
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &               IIPAR,     JJPAR,     1,         I0+1,
!     &               J0+1,      1,         REAL(ADJ) )
!
!      print*, 'MAKE_SAVE_FILE_3 : RAD   1 = ', CHK_RAD_1(IFD,JFD)
!      print*, 'MAKE_SAVE_FILE_3 : RAD   2 = ', CHK_RAD_2(IFD,JFD)
!      print*, 'MAKE_SAVE_FILE_3 : RAD   3 = ', CHK_RAD_3(IFD,JFD)
!      print*, 'MAKE_SAVE_FILE_3 : 2ord FD = ', TRACER_2D(IFD,JFD,1)
!      print*, 'MAKE_SAVE_FILE_3 : ADJ     = ', ADJ(IFD,JFD,1)
!
!      !========================================
!      ! Write ADJ / FD ratio
!      !========================================
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         IF ( ABS(TRACER_2D(I,J,1)) .gt. FILTER ) THEN 
!
!            TRACER_2D(I,J,1) = REAL(ADJ(I,J,1)) / TRACER_2D(I,J,1)
!          
!         ELSE 
!
!            TRACER_2D(I,J,1) = 1d0
!  
!         ENDIF 
!
!         ! Keep track of number of points that are off
!         IF ( ( TRACER_2D(I,J,1) > 1d0 + FD_DIFF ) .OR.
!     &        ( TRACER_2D(I,J,1) < 1D0 - FD_DIFF )      ) THEN 
!            N_OFF(I,J) = 1
!         ENDIF 
!
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      UNIT = 'none'
!
!         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &               HALFPOLAR, CENTER180, CATEGORY,  3,
!     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &               IIPAR,     JJPAR,     1,         I0+1,
!     &               J0+1,      1,         TRACER_2D )
!
!      ! print out statistics of the ADJ / FD ratio 
!      WRITE(6,*) '===================================================='
!      WRITE(6,*) ' Global validation test for values > ', FILTER
!      WRITE(6,*) ' MAX of global 2nd order ADJ / FD  = ',
!     &         MAXVAL(TRACER_2D(:,:,1))
!      WRITE(6,*) ' MIN of global 2nd order ADJ / FD  = ',
!     &         MINVAL(TRACER_2D(:,:,1))
!      WRITE(6,*) ' Number of places where ratio off by ',FD_DIFF,' = ',
!     &         SUM(N_OFF(:,:))
!      WRITE(6,*) '===================================================='
!
!      !========================================
!      ! Write 1st order FD gradient
!      !========================================
!      UNIT = 'kg/box'
!      !L = LFD
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         TRACER_2D(I,J,1) =(CHK_RAD_2(I,J) - CHK_RAD_1(I,J))
!     &          / FD_DIFF
!
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &            HALFPOLAR, CENTER180, CATEGORY,  4,
!     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &            IIPAR,     JJPAR,     1,         I0+1,
!     &            J0+1,      1,         TRACER_2D )
!
!      !========================================
!      ! Write chekpt values 
!      !========================================
!      UNIT = 'kg/box'
!      !L = LFD
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         TRACER_2D(I,J,1) = CHK_RAD_1(I,J)
!
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO
!
!      CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
!     &            HALFPOLAR, CENTER180, CATEGORY,  5,
!     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &            IIPAR,     JJPAR,     1,         I0+1,
!     &            J0+1,      1,         TRACER_2D )
!      ! Close file
!      CLOSE( IU_RST )
!      
!      !### Debug
!      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_SAVE_FILE_3: wrote file' )
!
!
!      ! Return to calling program
!      END SUBROUTINE MAKE_SAVE_FILE_3
!!----------------------------------------------------------------------

      SUBROUTINE MAKE_FD_FILE( YYYYMMDD, HHMMSS )
!
!******************************************************************************
!  Subroutine MAKE_FD_FILE creates GEOS-CHEM checkpt files of tracer 
!  concentrations [kg/box].
!  For use in checking chemistry adjoints.  (dkh, 07/19/06)  
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a checkpoint file       
!
!  NOTES:
! (1 ) Updated from MAKE_SAVE_FILE to v8, rename, replace CMN_ADJ, etc 
!       (dkh, ks, mak, cs, 06/09/09) 
!  
!******************************************************************************
!     
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : N_CALC 
      USE ADJ_ARRAYS_MOD,    ONLY : NFD
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJ_DIR
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,       ONLY : LPRT
      USE TIME_MOD,          ONLY : EXPAND_DATE, GET_TAU
      USE TRACER_MOD,        ONLY : N_TRACERS
      USE TRACER_MOD,        ONLY : STT

#     include "CMN_SIZE"          ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)        :: YYYYMMDD, HHMMSS

      ! Local Variables      
      INTEGER                    :: I,    I0, IOS, J,  J0, L, N
      INTEGER                    :: YYYY, MM, DD,  HH, SS
      CHARACTER(LEN=255)         :: FILENAME
      CHARACTER(LEN=255)         :: OUTPUT_FD_FILE

      ! Temporary storage arrays for checkpointed variables
      TYPE (XPLEX)                     :: TEMP(IIPAR,JJPAR,LLPAR)

      ! For binary punch file, version 2.0
      TYPE (XPLEX)                     :: LONRES, LATRES
      INTEGER, PARAMETER         :: HALFPOLAR = 1
      INTEGER, PARAMETER         :: CENTER180 = 1

      CHARACTER(LEN=20)          :: MODELNAME
      CHARACTER(LEN=40)          :: CATEGORY
      CHARACTER(LEN=40)          :: UNIT     
      CHARACTER(LEN=40)          :: RESERVED = ''
      CHARACTER(LEN=80)          :: TITLE 


      !=================================================================
      ! MAKE_FD_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      OUTPUT_FD_FILE = 'gctm.fd.YYYYMMDD.hhmm.NN'

      ! Clear some arrays 
      TEMP(:,:,:)  = 0e0

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM Checkpoint File: ' // 
     &           'Instantaneous Tracer Concentrations (v/v)'
      CATEGORY = 'IJ-CHK-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the save file for output -- binary punch format
      !=================================================================

      ! Copy the output checkpoint file name into a local variable
      FILENAME = TRIM( OUTPUT_FD_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Append the iteration number suffix to the file name
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add ADJ_DIR prefix to filename
      FILENAME = TRIM( DIAGADJ_DIR ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_FD_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      !=================================================================
      ! Write each checkpointed quantity to the checkpoint file
      !=================================================================
   
      ! Write the final concetration values as saved at the end of geos_mod.f
      UNIT = 'kg/box'

      ! Temporarily store data in CHECK_FINAL
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
            !TEMP(I,J,L) = CHK_STT(I,J,L,N)
         TEMP(I,J,L) = STT(I,J,L,NFD)
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  1,
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     LLPAR,     I0+1,
     &            J0+1,      1,         TEMP )

      ! Close file
      CLOSE( IU_RST )
      
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_FD_FILE: wrote file' )

      ! Return to calling program
      END SUBROUTINE MAKE_FD_FILE
!-----------------------------------------------------------------------

      SUBROUTINE MAKE_FDGLOB_FILE( YYYYMMDD, HHMMSS )
!
!******************************************************************************
!  Subroutine MAKE_FDGLOB_FILE creates GEOS-CHEM checkpt files of tracer 
!  concentrations [kg/box]. Like MAKE_FD_FILE, except calculate the finite
!  difference sensitivities directly.  Save these, the adjoint sensitivities, 
!  and the ratio adj / fd.  Save first and 2nd order finite difference 
!  sensitivities.  Requires running to XSTOP = 3. 
!  For use in checking process specific adjoints. (dkh, 01/23/07) 
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a checkpoint file       
!
!  NOTES:
!  (1 ) Now write out both sets of 1st order FD gradients, and define a new 
!        category (FD-TEST) for viewing in gamap.    (dkh, 10/10/08) 
!  (2 ) Updated from MAKE_SAVE_FILE to v8, rename, replace CMN_ADJ, etc 
!       (dkh, ks, mak, cs, 06/09/09) 
!  (3 ) NNEMS Now in tracerid_adj_mod.f (mak, 6/14/09)
!  (4 ) Updated to include LADJ_STRAT (hml, dkh, 02/14/12, adj32_025) 
!******************************************************************************
!     
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,    ONLY : NFD, EMSFD, MFD, LFD, MMSCL
      USE ADJ_ARRAYS_MOD,    ONLY : FD_DIFF, NNEMS
      USE ADJ_ARRAYS_MOD,    ONLY : ICSFD
      USE ADJ_ARRAYS_MOD,    ONLY : IFD, JFD
      USE ADJ_ARRAYS_MOD,    ONLY : STRFD
      USE ADJ_ARRAYS_MOD,    ONLY : NSTPL
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJ_DIR
      USE DIRECTORY_ADJ_MOD, ONLY : OPTDATA_DIR
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR
      USE ERROR_MOD,         ONLY : DEBUG_MSG,   ERROR_STOP
      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,       ONLY : LPRT 
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_EMS, LICS
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_STRAT
      USE TIME_MOD,          ONLY : EXPAND_DATE, GET_TAU
      USE TRACER_MOD,        ONLY : N_TRACERS

#     include "CMN_SIZE"          ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)        :: YYYYMMDD, HHMMSS

      ! Local Variables      
      INTEGER                    :: I, I0, IOS, J,  J0, L, N
      INTEGER                    :: YYYY, MM, DD,  HH, SS
      CHARACTER(LEN=255)         :: FILENAME
      CHARACTER(LEN=255)         :: INPUT_GDT_FILE
      CHARACTER(LEN=255)         :: INPUT_FD_FILE
      CHARACTER(LEN=255)         :: OUTPUT_FDGLOB_FILE
      INTEGER                    :: N_OFF(IIPAR,JJPAR)
      TYPE (XPLEX), PARAMETER          :: FILTER = xplex(1d0,0d0)
 
      ! Temporary storage arrays for checkpointed variables
      TYPE (XPLEX)                     :: TEMP1(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                     :: TEMP2(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                     :: TEMP3(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                     :: TRACER(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                     :: EMS_3D(IIPAR,JJPAR,MMSCL)
      TYPE (XPLEX)                     :: ADJ_2D(IIPAR,JJPAR,1)
      TYPE (XPLEX)                     :: TRACER_2D(IIPAR,JJPAR,1)
      COMPLEX*16                     :: D_TEMP1(IIPAR,JJPAR,LLPAR)
      COMPLEX*16                     :: D_TEMP2(IIPAR,JJPAR,LLPAR)
      COMPLEX*16                     :: D_TEMP3(IIPAR,JJPAR,LLPAR)
      COMPLEX*16                     :: D_TRACER(IIPAR,JJPAR,LLPAR)
      COMPLEX*16                     :: D_EMS_3D(IIPAR,JJPAR,MMSCL)
      COMPLEX*16                     :: D_ADJ_2D(IIPAR,JJPAR,1)
      COMPLEX*16                     :: D_TRACER_2D(IIPAR,JJPAR,1)

      ! For strat prod and loss (hml, 09/01/11, adj32_025)
      TYPE (XPLEX)                     :: PROD_3D(IIPAR,JJPAR,MMSCL)
      TYPE (XPLEX)                     :: LOSS_3D(IIPAR,JJPAR,MMSCL)
      TYPE (XPLEX)                     :: EMS_N_3D(IIPAR,JJPAR,MMSCL)
      COMPLEX*16                     :: D_PROD_3D(IIPAR,JJPAR,MMSCL)
      COMPLEX*16                     :: D_LOSS_3D(IIPAR,JJPAR,MMSCL)
      COMPLEX*16                     :: D_EMS_N_3D(IIPAR,JJPAR,MMSCL)
      ! For binary punch file, version 2.0
      TYPE (XPLEX)                     :: LONRES, LATRES
      COMPLEX*16                     :: D_LONRES, D_LATRES
      INTEGER                    :: HALFPOLAR
      INTEGER                    :: CENTER180

      CHARACTER(LEN=20)          :: MODELNAME
      CHARACTER(LEN=40)          :: CATEGORY
      CHARACTER(LEN=40)          :: UNIT     
      CHARACTER(LEN=40)          :: RESERVED = ''
      CHARACTER(LEN=80)          :: TITLE 

      INTEGER                    :: NI,     NJ,     NL
      INTEGER                    :: IFIRST, JFIRST, LFIRST
      INTEGER                    :: NTRACER,   NSKIP
      TYPE (XPLEX)                     :: ZTAU0,     ZTAU1
      COMPLEX*16                     :: D_ZTAU0,     D_ZTAU1

      !=================================================================
      ! MAKE_FDGLOB_FILE begins here!
      !=================================================================

      ! Clear some arrays 
      TEMP1(:,:,:)        = 0e0
      TEMP2(:,:,:)        = 0e0
      TEMP3(:,:,:)        = 0e0
      EMS_3D(:,:,:)       = 0d0
      ADJ_2D(:,:,:)       = 0e0
      TRACER(:,:,:)       = 0e0
      TRACER_2D(:,:,:)    = 0e0
      N_OFF(:,:)          = 0d0
      D_TEMP1(:,:,:)        = 0e0
      D_TEMP2(:,:,:)        = 0e0
      D_TEMP3(:,:,:)        = 0e0
      D_EMS_3D(:,:,:)       = 0d0
      D_ADJ_2D(:,:,:)       = 0e0
      D_TRACER(:,:,:)       = 0e0
      D_TRACER_2D(:,:,:)    = 0e0
      
      ! strat prod and loss (hml)
      PROD_3D(:,:,:)      = 0d0
      LOSS_3D(:,:,:)      = 0d0
      EMS_N_3D(:,:,:)     = 0d0
      D_PROD_3D(:,:,:)      = 0d0
      D_LOSS_3D(:,:,:)      = 0d0
      D_EMS_N_3D(:,:,:)     = 0d0

      !========================================
      ! Read *.fd.* file from unperturbed run
      !========================================
      INPUT_FD_FILE = 'gctm.fd.YYYYMMDD.hhmm.NN'

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_FD_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Append the iteration number suffix to the file name
      CALL EXPAND_NAME( FILENAME, 1 )

      ! Add ADJ_DIR prefix to name
      FILENAME = TRIM( DIAGADJ_DIR ) // TRIM( FILENAME )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'F D    F I L E   I N P U T: base '


      WRITE( 6, 400 ) TRIM( FILENAME )
 400  FORMAT( '     - READ_FDGLOB_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
      
      !=================================================================
      ! Read checkpointed variables 
      !=================================================================
      READ( IU_RST, IOSTAT=IOS )
     &  MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
       LATRES%r = dble(D_LATRES)
       LONRES%r = dble(D_LONRES)
       LATRES%i = dimag(D_LATRES)
       LONRES%i = dimag(D_LONRES)
      !! IOS < 0 is end-of-file, so exit
      !IF ( IOS < 0 ) EXIT

      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 ) 
     &   CALL IOERROR( IOS,IU_RST,'read_fdglob_file:1' )

      READ( IU_RST, IOSTAT=IOS )
     &      CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &      NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &      NSKIP
       ZTAU0%r = dble(D_ZTAU0)
       ZTAU1%r = dble(D_ZTAU1)
       ZTAU0%i = dimag(D_ZTAU0)
       ZTAU1%i = dimag(D_ZTAU1)
      IF ( IOS /= 0 ) 
     &   CALL IOERROR(IOS,IU_RST,'read_fdglob_file:2' )

      READ( IU_RST, IOSTAT=IOS )
     &     ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
       TRACER%r = dble(D_TRACER)
       TRACER%i = dimag(D_TRACER)
      IF ( IOS /= 0 ) 
     &   CALL IOERROR( IOS,IU_RST,'read_fdglob_file:3' )

      !==============================================================
      ! Assign data from the TRACER array to the STT array.
      !==============================================================

      ! Only process checkpoint data (i.e. mixing ratio)
      IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN

          ! Make sure array dimensions are of global size
          ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
          !CALL CHECK_DIMENSIONS( NI, NJ, NL )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
          DO L = 1, LLPAR
          DO J = 1, JJPAR
          DO I = 1, IIPAR
             TEMP1(I,J,L) = TRACER(I,J,L)
          ENDDO
          ENDDO
          ENDDO
!$OMP END PARALLEL DO

      ENDIF

      ! Close file
      CLOSE( IU_RST )
      

      !========================================
      ! Read *.fd* file from perturbed run
      !========================================
      INPUT_FD_FILE = 'gctm.fd.YYYYMMDD.hhmm.NN'

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_FD_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Append the iteration number suffix to the file name
      CALL EXPAND_NAME( FILENAME, 2 )

      ! Add ADJ_DIR prefix to name
      FILENAME = TRIM( DIAGADJ_DIR ) // TRIM( FILENAME )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'F D    F I L E   I N P U T: +pert '


      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_FDGLOB_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
      
      !=================================================================
      ! Read checkpointed variables 
      !=================================================================
      ! Read the values of STT
      READ( IU_RST, IOSTAT=IOS )
     &  MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
        LONRES%r = dble(D_LONRES)
        LATRES%r = dble(D_LATRES)
        LONRES%i = dimag(D_LONRES)
        LATRES%i = dimag(D_LATRES)
      !! IOS < 0 is end-of-file, so exit
      !IF ( IOS < 0 ) EXIT

      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 ) 
     &   CALL IOERROR( IOS,IU_RST,'read_fdglob_file:4' )

      READ( IU_RST, IOSTAT=IOS )
     &      CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &      NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &      NSKIP
       ZTAU0%r = dble(D_ZTAU0)
       ZTAU1%r = dble(D_ZTAU1)
       ZTAU0%i = dimag(D_ZTAU0)
       ZTAU1%i = dimag(D_ZTAU1)
      IF ( IOS /= 0 ) 
     &   CALL IOERROR(IOS,IU_RST,'read_fdglob_file:5' )

      READ( IU_RST, IOSTAT=IOS )
     &     ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
        TRACER%r = dble(D_TRACER)
        TRACER%i = dimag(D_TRACER)
      IF ( IOS /= 0 ) 
     &    CALL IOERROR( IOS,IU_RST,'read_fdglob_file:6' )

      !==============================================================
      ! Assign data from the TRACER array to the STT array.
      !==============================================================

      ! Only process checkpoint data (i.e. mixing ratio)
      IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN

         ! Make sure array dimensions are of global size
         ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
         !CALL CHECK_DIMENSIONS( NI, NJ, NL )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            TEMP2(I,J,L) = TRACER(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF
 
      ! Close file
      CLOSE( IU_RST )

      !========================================
      ! Read *.fd.* file from 2nd perturbed run
      !========================================
      INPUT_FD_FILE = 'gctm.fd.YYYYMMDD.hhmm.NN'

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_FD_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Append the iteration number suffix to the file name
      CALL EXPAND_NAME( FILENAME, 3 )

      ! Add ADJ_DIR prefix to name
      FILENAME = TRIM( DIAGADJ_DIR ) // TRIM( FILENAME )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'F D    F I L E   I N P U T: -pert '


      WRITE( 6, 888 ) TRIM( FILENAME )
 888  FORMAT( '     - READ_FDGLOB_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
      
      !=================================================================
      ! Read checkpointed variables 
      !=================================================================
      ! Read the values of CHK_STT
      READ( IU_RST, IOSTAT=IOS )
     &  MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
       LONRES%r = dble(D_LONRES)
        LATRES%r = dble(D_LATRES)
        LONRES%i = dimag(D_LONRES)
        LATRES%i = dimag(D_LATRES)
      !! IOS < 0 is end-of-file, so exit
      !IF ( IOS < 0 ) EXIT

      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 ) 
     &   CALL IOERROR( IOS,IU_RST,'read_fdglob_file:7' )

      READ( IU_RST, IOSTAT=IOS )
     &      CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &      NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &      NSKIP
        ZTAU0%r = dble(D_ZTAU0)
        ZTAU1%r = dble(D_ZTAU1)
        ZTAU0%i = dimag(D_ZTAU0)
        ZTAU1%i = dimag(D_ZTAU1)
      IF ( IOS /= 0 ) 
     &   CALL IOERROR(IOS,IU_RST,'read_fdglob_file:8' )

      READ( IU_RST, IOSTAT=IOS )
     &     ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
        TRACER%r = dble(D_TRACER)
        TRACER%i = dimag(D_TRACER)
      IF ( IOS /= 0 ) 
     &   CALL IOERROR( IOS,IU_RST,'read_fdglob_file:9' )

      !==============================================================
      ! Assign data from the TRACER array to the STT array.
      !==============================================================

      ! Only process checkpoint data (i.e. mixing ratio)
      IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN

         ! Make sure array dimensions are of global size
         ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
         !CALL CHECK_DIMENSIONS( NI, NJ, NL )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            TEMP3(I,J,L) = TRACER(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF
      
      ! Close file
      CLOSE( IU_RST )

      !========================================
      ! Read *.gdt.* file from unperturbed run
      !========================================
      INPUT_GDT_FILE = 'gctm.gdt.01'

      ! Initialize some variables
      TRACER(:,:,:) = 0e0

      !=================================================================
      ! Open gradient file and read top-of-file header
      !=================================================================

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_GDT_FILE )

      ! Add OPT_DATA_DIR prefix to FILENAME
      FILENAME = TRIM( OPTDATA_DIR ) // TRIM( FILENAME )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'G D T   F I L E   I N P U T'
      WRITE( 6, 101 ) TRIM( FILENAME )
 101  FORMAT( 'READ_GDT_FILE: Reading ', a )

      
      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )

      !=================================================================
      ! Read adjoints -- store in the TRACER array
      !=================================================================

      IF ( LADJ_EMS ) THEN 
         DO N = 1, NNEMS
            READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) EXIT

            ! IOS > 0 is a real I/O error -- print error message
            IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'fdglob_file:10' )

            READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP

            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'fdglob_file:11')

            READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( EMS_3D(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'fdglob_file:12')

            ! Don't write if LADJ_STRAT (hml, 10/31/11, adj32_025)
            IF ( .NOT. LADJ_STRAT ) THEN

               !==============================================================
               ! Assign data from the TRACER array to the ADJ_STT array.
               !==============================================================

               ! Save the gradients selected by EMSFD and MFD 
               IF ( CATEGORY(1:8) == 'IJ-GDE-$' .and. N == EMSFD ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
                  DO J = 1, JJPAR
                  DO I = 1, IIPAR
                     ADJ_2D(I,J,1) = EMS_3D(I,J,MFD)
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO

               ENDIF
            ENDIF
         ENDDO

         ! Read GDEN (hml, 09/11/11, adj32_025)
         DO N = 1, NNEMS

            READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) EXIT

            ! IOS > 0 is a real I/O error -- print error message
            IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'fdglob_file:10-b' )

            READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP

            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'fdglob_file:11-b')

            READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( EMS_N_3D(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'fdglob_file:12-b')

         ENDDO

         ! For strat prod and loss (hml, 08/30/11, adj32_025)
         IF ( LADJ_STRAT ) THEN

            ! Strat production
            DO N = 1, NSTPL

               READ( IU_RST, IOSTAT=IOS )
     &           MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

               ! IOS < 0 is end-of-file, so exit
               IF ( IOS < 0 ) EXIT

               ! IOS > 0 is a real I/O error -- print error message
               IF ( IOS > 0 ) CALL IOERROR(IOS,IU_RST,'fdglob_file:13')

               READ( IU_RST, IOSTAT=IOS )
     &              CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &              NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &              NSKIP

               IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_RST,'fdglob_file:14')

               READ( IU_RST, IOSTAT=IOS )
     &              ( ( ( PROD_3D(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

               IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_RST,'fdglob_file:15')

            ENDDO

            ! Strat loss
            DO N = 1, NSTPL

               READ( IU_RST, IOSTAT=IOS )
     &           MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

               ! IOS < 0 is end-of-file, so exit
               IF ( IOS < 0 ) EXIT

               ! IOS > 0 is a real I/O error -- print error message
               IF ( IOS > 0 ) CALL IOERROR(IOS,IU_RST,'fdglob_file:16')

               READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP

               IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_RST,'fdglob_file:17')

               READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( LOSS_3D(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

               IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_RST,'fdglob_file:18')

               !========================================================
               ! Assign data from the LOSS_3D array to the ADJ_STT array.
               !========================================================

               ! Save the gradients selected by EMSFD and MFD 
               IF ( CATEGORY(1:8) == 'IJ-GDL-$' .and. N == STRFD ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
                  DO J = 1, JJPAR
                  DO I = 1, IIPAR
                     ADJ_2D(I,J,1) = LOSS_3D(I,J,MFD)
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO

               ENDIF
            ENDDO
         ENDIF

      ELSEIF ( LICS ) THEN 

         TRACER(:,:,:) = 0d0 

         DO N = 1, N_TRACERS
            READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) EXIT

            ! IOS > 0 is a real I/O error -- print error message
            IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'fdglob_file:13' )

            READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP

            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'fdglob_file:14')

            READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'fdglob_file:15')

            !==============================================================
            ! Assign data from the TRACER array to the ADJ_STT array.
            !==============================================================

            ! Save the gradients selected by EMSFD and MFD 
            IF ( CATEGORY(1:8) == 'IJ-GDT-$' .and. N == ICSFD ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  ADJ_2D(I,J,1) = TRACER(I,J,LFD)
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

            ENDIF
         ENDDO

      ENDIF 

      ! Close file
      CLOSE( IU_RST )


      !========================================
      ! Write 2nd order FD gradient
      !========================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         TRACER_2D(I,J,1) =(TEMP2(I,J,LFD) - TEMP3(I,J,LFD))
     &          /  ( 2d0 * FD_DIFF )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Hardwire output file for now
      OUTPUT_FDGLOB_FILE = 'gctm.fdglob.YYYYMMDD.hhmm'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM Checkpoint File: ' // 
     &           'Instantaneous Tracer Concentrations (v/v)'
      CATEGORY = 'FD-TEST'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the save file for output -- binary punch format
      !=================================================================

      ! Copy the output checkpoint file name into a local variable
      FILENAME = TRIM( OUTPUT_FDGLOB_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to filename
      FILENAME = TRIM( DIAGADJ_DIR ) // TRIM( FILENAME )

      WRITE( 6, 102 ) TRIM( FILENAME )
 102  FORMAT( '     - MAKE_FDGLOB_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      UNIT = 'kg/box'

         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  1,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     1,         I0+1,
     &               J0+1,      1,         TRACER_2D )

      !========================================
      ! Write ADJ gradient 
      !========================================
      UNIT = 'kg/box'

         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  2,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     1,         I0+1,
     &               J0+1,      1,         ADJ_2D )

      !========================================
      ! Write ADJ / FD ratio
      !========================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         IF ( ABS(TRACER_2D(I,J,1)) .gt. FILTER ) THEN 

            TRACER_2D(I,J,1) = ADJ_2D(I,J,1) / TRACER_2D(I,J,1)
          
         ELSE 

            TRACER_2D(I,J,1) = 1d0
  
         ENDIF 

         ! Keep track of number of points that are off
         IF ( ( TRACER_2D(I,J,1) > 1d0 + FD_DIFF ) .OR.
     &        ( TRACER_2D(I,J,1) < 1D0 - FD_DIFF )      ) THEN 
            N_OFF(I,J) = 1
         ENDIF 

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      UNIT = 'none'

         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  3,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     1,         I0+1,
     &               J0+1,      1,         TRACER_2D )

      ! print out statistics of the ADJ / FD ratio 
      WRITE(6,*) '===================================================='
      WRITE(6,*) ' Global validation test for values > ', FILTER
      WRITE(6,*) ' MAX of global 2nd order ADJ / FD  = ',
     &         MAXVAL(TRACER_2D(:,:,1)), MAXLOC((TRACER_2D(:,:,1)%r))
      WRITE(6,*) ' MIN of global 2nd order ADJ / FD  = ',
     &         MINVAL(TRACER_2D(:,:,1)), MINLOC((TRACER_2D(:,:,1)%r))
      WRITE(6,*) ' Number of places where ratio off by ',FD_DIFF,' = ',
     &         SUM(N_OFF(:,:))
      WRITE(6,*) '===================================================='

      print*, ' FD2 ' , (TEMP2(IFD,JFD,LFD) - TEMP3(IFD,JFD,LFD))
     &          /  ( 2d0 * FD_DIFF )
      print*, ' FD2 from ', TEMP2(IFD,JFD,LFD) , TEMP3(IFD,JFD,LFD)
      print*, ' ADJ ' , ADJ_2D(IFD,JFD,1)

      !========================================
      ! Write 1st order FD gradient
      !========================================
      UNIT = 'kg/box'
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         TRACER_2D(I,J,1) =( TEMP2(I,J,LFD) - TEMP1(I,J,LFD) )
     &          / FD_DIFF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  4,
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     1,         I0+1,
     &            J0+1,      1,         TRACER_2D )

      !========================================
      ! Write chekpt values 
      !========================================
      UNIT = 'kg/box'
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         TRACER(I,J,1) = TEMP1(I,J,LFD)
         TRACER(I,J,2) = TEMP2(I,J,LFD)
         TRACER(I,J,3) = TEMP3(I,J,LFD)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  5,
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     3,         I0+1,
     &            J0+1,      1,     (TRACER(:,:,1:3)) )

      !========================================
      ! Write the other 1st order FD gradient. (dkh, 10/10/08) 
      !========================================
      UNIT = 'kg/box'
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         TRACER_2D(I,J,1) =( TEMP3(I,J,LFD) - TEMP1(I,J,LFD) )
     &          / ( - FD_DIFF )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  6,
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     1,         I0+1,
     &            J0+1,      1,         TRACER_2D )

      ! Close file
      CLOSE( IU_RST )
      
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_FDGLOB_FILE: wrote file' )


      ! Return to calling program
      END SUBROUTINE MAKE_FDGLOB_FILE
!----------------------------------------------------------------------

      SUBROUTINE MAKE_CHK_CON_FILE( YYYYMMDD, HHMMSS )
!
!******************************************************************************
!  Subroutine MAKE_CHK_CON_FILE creates GEOS-CHEM checkpt files of tracer 
!  mixing ratios (v/v), and exit values in binary punch file format. 
!  (mak, 8/2/07)
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a checkpoint file       
!
!  Passed via CMN:
!  ============================================================================
!  (1 ) TAU    : TAU value (elapsed hours) at start of diagnostic interval
!
!  Passed via ???:
!  ============================================================================
!  (1 ) CHECKPT	: Array of quantities to be checkpointed     
!                               dim=(IIPAR,JJPAR,LLPAR,NCHECKPT) 
!
!  NOTES:
!   Just like MAKE_RESTART_FILE except
!	- only include quantities used as input to RPMARES
!       - include hhmmss in file name
!       - writes files to ADJ_DIR and can zip them
!       dkh, 9/30/04
!  (2 ) Zip *.chk.* files one day at a time in a parallel loop. Add access
!        to GET_TS_CHEM.   (dkh, 11/22/04)  
!  (3 ) Add support for L_RECOMP option to recompute (rather than checkpoint)
!        variables RP_OUT etc.  (dkh, 02/09/05)
!  (4 ) Now write values from CHK_STT_CON. (mak, 8/2/07)
!  (5 ) Change file names to *.chk.con.* so they get cleaned out by shell scripts
!        that purge *.chk.* files. (dkh, 10/10/08) 
!  (6 ) Update to v8 adj (dkh, 06/11/09) 
!******************************************************************************
!     
      ! References to F90 modules
      USE BPCH2_MOD
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP, ALLOC_ERR
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR
      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,       ONLY : LPRT 
      USE TIME_MOD,          ONLY : EXPAND_DATE, GET_TAU
      USE TRACER_MOD,        ONLY : ITS_A_TAGCO_SIM
      USE TRACER_MOD,        ONLY : N_TRACERS

#     include "CMN_SIZE"          ! Size parameters
#     include "comode.h"          ! IGAS

      ! Arguments
      INTEGER, INTENT(IN)        :: YYYYMMDD, HHMMSS

      ! Local Variables      
      INTEGER              :: I,    I0, IOS, J,  J0, L, N, JLOOP
      INTEGER              :: YYYY, MM, DD,  HH, SS, ZIP_HH, AS
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      TYPE (XPLEX)               :: LONRES, LATRES
      TYPE (XPLEX)               :: CHECK_FINAL(IIPAR,JJPAR,LLPAR)
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      INTEGER 		   :: MAX_nitr_max
      INTEGER 		   :: NSOFAR

      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 
      TYPE (XPLEX)               :: nitr_max_real(IIPAR, JJPAR, LLPAR)

      !=================================================================
      ! MAKE_CHK_CON_FILE begins here!
      !=================================================================

      ! NEW: rename them *.chk.con.*  (dkh, 10/10/08) 
      OUTPUT_CHECKPT_FILE = 'gctm.chk.con.YYYYMMDD.hhmm'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM Convection Checkpoint File: ' // 
     &           'Instantaneous Tracer Concentrations (v/v)'
      CATEGORY = 'IJ-CHK-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the checkpoint file for output -- binary punch format
      !=================================================================

      ! Copy the output checkpoint file name into a local variable
      FILENAME = TRIM( OUTPUT_CHECKPT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to filename
      FILENAME = TRIM( ADJTMP_DIR ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_CHECKPT_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
      !=================================================================
      ! Write each checkpointed quantity to the checkpoint file
      !=================================================================
    
      IF ( ITS_A_TAGCO_SIM() )THEN
         UNIT = 'v/v'
      ENDIF

      DO N = 1, N_TRACERS

         ! Temporarily store data in CHECK_FINAL
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            CHECK_FINAL(I,J,L) = CHK_STT_CON(I,J,L,N)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     LLPAR,     I0+1,
     &               J0+1,      1,         CHECK_FINAL )
      ENDDO


      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_CHK_CON_FILE: wrote file' )

      ! Return to calling program
      END SUBROUTINE MAKE_CHK_CON_FILE

!------------------------------------------------------------------------------
      SUBROUTINE READ_CHK_CON_FILE( YYYYMMDD, HHMMSS ) 
!
!******************************************************************************
!  Subroutine READ_CHK_CON_FILE initializes GEOS-CHEM tracer concentrations 
!  from a checkpoint file (binary punch file format) from before convection
!  (dkh, 8/30/04, mak, 8/2/07)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!
!  Passed via CMN:
!  ============================================================================
!  (1 ) TAU    : TAU value (elapsed hours) at start of diagnostic interval
!
!  Notes
!  (1 ) Just like READ_RESTART_FILE except
!	- load the variables from TRACER directly back into the CHECKPT array
!       - file name now includes hhmmss
!       - reads files from ADJ_DIR (and can unzip them if L_ZIP_CHECKPT) 
!       - removes .chk. files after reading (if L_DEL_CHECKPT)
!       dkh, 9/30/04
!  (2 ) Add DATE(2) and reference GET_NHMDe and GET_NHMSe to enable BATCH_ZIP
!        (dkh, 11/22/04)
!  (3 ) Read in CHK_STT_CON (mak, 8/2/07)
!  (4 ) Rename from *.chkcon.* to *.chk.con.* (dkh, 10/10/08) 
!  (5 ) Delete files after they've been used. (dkh, 10/10/08) 
!  (6 ) Remove the IF ( N == 1 ) line. (dkh, 10/10/08) 
!  (7 ) Remove obsolete options (L_DEL_CHECKPT, L_ZIP_CHECKPT, L_RECOMP), 
!        check for aeroosl simulation (LSULF) and update names to v8 (dkh, 06/11/09) 
!  (8 ) Update to v8 adj (dkh, 06/11/09) 
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,         ONLY : OPEN_BPCH2_FOR_READ
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR 
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP, ALLOC_ERR
      USE FILE_MOD,          ONLY : IU_RST, IOERROR
      USE LOGICAL_MOD,       ONLY : LPRT
      USE LOGICAL_ADJ_MOD,   ONLY : LDEL_CHKPT
      USE RESTART_MOD,       ONLY : CHECK_DIMENSIONS
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TRACER_MOD,        ONLY : N_TRACERS
      USE UNIX_CMDS_MOD,     ONLY : REMOVE_CMD


#     include "CMN_SIZE"   ! Size parameters
#     include "comode.h"   ! ITLOOP, IGAS

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER             :: I, IOS, J, L, N, JLOOP, NN, NTL, AS
      INTEGER             :: NCOUNT(NNPAR) 
      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
      COMPLEX*16              :: D_TRACER(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME
      CHARACTER(LEN=255)  :: REMOVE_CHK_FILE_CMD
     

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      TYPE (XPLEX)              :: LONRES,    LATRES
      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
      COMPLEX*16              :: D_LONRES,    D_LATRES
      COMPLEX*16              :: D_ZTAU0,     D_ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED

      !=================================================================
      ! READ_CHECKPT_FILE begins here!
      !=================================================================

      INPUT_CHECKPT_FILE = 'gctm.chk.con.YYYYMMDD.hhmm'

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0
      D_TRACER(:,:,:) = 0e0
      !=================================================================
      ! Open checkpoint file and read top-of-file header
      !=================================================================
      
      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_CHECKPT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to name
      FILENAME = TRIM( ADJTMP_DIR ) // TRIM( FILENAME )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'C H E C K P T   F I L E   I N P U T'

      ! Remove obsolete option
!      ! Unzip checkpt file
!      IF ( L_ZIP_CHECKPT ) THEN
!         UNZIP_FILE_CMD = TRIM( GUNZIP_CMD ) // ' ' //    
!     &                  TRIM( FILENAME ) // ZIP_SUFFIX    
!         CALL SYSTEM( TRIM( UNZIP_FILE_CMD ) )
!         WRITE( 6, 99 ) TRIM( UNZIP_FILE_CMD ) 
! 99   FORMAT( '     - READ_CHECKPT_FILE: Executing: ',a )
!      ENDIF


      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_CHECKPT_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
      
      !=================================================================
      ! Read checkpointed variables 
      !=================================================================
 
      ! Read the values of CHK_STT_CON
      !DO N = 1, NOBS!+1
      DO N = 1, N_TRACERS
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LONRES%r = dble(D_LONRES)
         LATRES%r = dble(D_LATRES)
         LONRES%i = dimag(D_LONRES)
         LATRES%i = dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:7' )

         READ( IU_RST, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0%r = dble(D_ZTAU0)
         ZTAU1%r = dble(D_ZTAU1)
         ZTAU0%i = dimag(D_ZTAU0)
         ZTAU1%i = dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:8' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
         TRACER%r = dble(D_TRACER)
         TRACER%i = dimag(D_TRACER)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:9' )

         !==============================================================
         ! Assign data from the TRACER array to the STT array.
         !==============================================================

         ! Only process checkpoint data (i.e. mixing ratio)
         IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN

            ! Make sure array dimensions are of global size
            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
            !print*, 'before check_dimensions ni, nj, nl are', ni, nj, nl
            CALL CHECK_DIMENSIONS( NI, NJ, NL )

               
            ! Remove (dkh, 10/10/08) 
            ! IF ( N == 1) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  CHK_STT_CON(I,J,L,N) = TRACER(I,J,L)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

            !ENDIF 

         ENDIF !category is checkpoint
      ENDDO


      ! Close file
      CLOSE( IU_RST )      


 555  CONTINUE

      ! Remove files if L_CHK_DEL = TRUE 
      IF ( LDEL_CHKPT ) THEN 
         REMOVE_CHK_FILE_CMD  = TRIM ( REMOVE_CMD ) // ' ' //
     &                          TRIM ( FILENAME )

        CALL SYSTEM( TRIM( REMOVE_CHK_FILE_CMD ) )

        WRITE( 6, 102 ) TRIM( REMOVE_CHK_FILE_CMD )
 102    FORMAT( '     - READ_CHECKPT_FILE: Executing: ',a )
      ENDIF 


      ! Remove obsolete  (dkh, 06/11/09) 
!      ! Zip the .chk. file if it hasn't been deleted and zipping 
!      ! is requested
!      IF ( L_ZIP_CHECKPT .AND. (.NOT. L_DEL_CHECKPT) ) THEN
!         CALL BATCH_ZIP( YYYYMMDD, HHMMSS, 'chk', -1 )
!      ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_CHK_CON_FILE: read file' )

      ! Return to calling program
      END SUBROUTINE READ_CHK_CON_FILE

!-----------------------------------------------------------------------

      SUBROUTINE MAKE_CHK_DYN_FILE( YYYYMMDD, HHMMSS )
!
!******************************************************************************
!  Subroutine MAKE_CHK_DYN_FILE creates GEOS-CHEM checkpt files 
!  at the dynamic time step. (dkh, 02/01/09) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create a checkpoint file       
!
!  Passed via CMN:
!  ============================================================================
!  (1 ) TAU    : TAU value (elapsed hours) at start of diagnostic interval
!
!  Passed via module variables
!  ============================================================================
!  (1 ) to add....
!                
!
!  NOTES:
!   Just like MAKE_CHK_FILE except
!
!  (1 ) Now checkpoint T_DAY, T_15_AVG  (dkh, 01/23/10)
!******************************************************************************
!     
      ! References to F90 modules
      USE BPCH2_MOD,         ONLY : BPCH3,     GET_MODELNAME
      USE BPCH2_MOD,         ONLY : BPCH2,     OPEN_BPCH2_FOR_WRITE
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR 
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP, ALLOC_ERR
      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE MEGAN_MOD,         ONLY : GET_T_DAY
      USE MEGAN_MOD,         ONLY : GET_T_15_AVG
      USE MEGAN_MOD,         ONLY : DAY_DIM
      USE MEGAN_MOD,         ONLY : CHK_T_15_AVG
      USE MEGAN_MOD,         ONLY : CHK_T_DAY
      USE TIME_MOD,          ONLY : EXPAND_DATE, GET_TAU
      USE TIME_MOD,          ONLY : GET_TS_CONV
      USE TIME_MOD,          ONLY : ITS_TIME_FOR_A3_ADJ
      USE TIME_MOD,          ONLY : ITS_TIME_TO_CHK_T_15_AVG
      USE LOGICAL_MOD,       ONLY : LWETD, LCONV
      USE LOGICAL_MOD,       ONLY : LPRT
      USE LOGICAL_MOD,       ONLY : LMEGAN 
      USE LOGICAL_ADJ_MOD,   ONLY : LAERO_THERM
      USE TRACER_MOD,        ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,        ONLY : ITS_AN_AEROSOL_SIM


#     include "CMN_SIZE"   ! Size parameters
#     include "comode.h"   ! IGAS

      ! Arguments
      INTEGER, INTENT(IN)  :: YYYYMMDD, HHMMSS

      ! Local Variables      
      INTEGER              :: I,    I0, IOS, J,  J0, L, N, JLOOP
      INTEGER              :: YYYY, MM, DD,  HH, SS, ZIP_HH, AS
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      TYPE (XPLEX)               :: LONRES, LATRES
      TYPE (XPLEX)               :: CHECK_FINAL(IIPAR,JJPAR,LLPAR)
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      INTEGER 		   :: MAX_nitr_max
      INTEGER 		   :: NSOFAR
      INTEGER 		   :: NS
      INTEGER 		   :: NSTEP
      INTEGER 		   :: CONVDT

      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 
      TYPE (XPLEX)               :: nitr_max_real(IIPAR, JJPAR, LLPAR)

      !=================================================================
      ! MAKE_CHK_DYN_FILE begins here!
      !=================================================================
      OUTPUT_CHECKPT_FILE = 'gctm.chk.dyn.YYYYMMDD.hhmm'


      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM Convection Checkpoint File' 
      CATEGORY = 'IJ-CHKD$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the checkpoint file for output -- binary punch format
      !=================================================================

      ! Copy the output checkpoint file name into a local variable
      FILENAME = TRIM( OUTPUT_CHECKPT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to filename
      FILENAME = TRIM( ADJTMP_DIR ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_CHK_DYN_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      ! Initialize NSOFAR
      NSOFAR = 0 
      
      !=================================================================
      ! Write each checkpointed quantity to the checkpoint file
      !=================================================================

      Unit = 'hPa'

      ! Write the surface pressures before and after transport (IIPAR,JJPAR,2)
      CALL BPCH3( IU_RST,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  NSOFAR + 1,
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     2,         I0+1,
     &            J0+1,      1,         CHK_PSC)

      ! Set NSOFAR
      NSOFAR = NSOFAR + 1


      IF ( LWETD .and. 
     &   ( ITS_AN_AEROSOL_SIM() .or. ITS_A_FULLCHEM_SIM() ) ) THEN



         ! H2O2s 
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  1 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     LLTROP,     I0+1,
     &               J0+1,      1,         WETD_CHK_H2O2s )

         ! SO2s
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  2 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     LLTROP,     I0+1,
     &               J0+1,      1,         WETD_CHK_SO2s )

         ! SO4
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  3 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     LLTROP,     I0+1,
     &               J0+1,      1,         WETD_CHK_SO4 )

         ! SO2 
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  4 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     LLTROP,     I0+1,
     &               J0+1,      1,         WETD_CHK_SO2 )

         NSOFAR = NSOFAR + 4 

      ENDIF  ! LWETD


      ! Write the concentrations used in convection
      IF ( LCONV .AND. 
     &   ( ITS_AN_AEROSOL_SIM() .or. ITS_A_FULLCHEM_SIM() ) ) THEN

         UNIT = 'kg'

         ! H2O2s 
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  1 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     LLTROP,     I0+1,
     &               J0+1,      1,         CONV_CHK_H2O2s )

         ! SO2s 
         CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  2 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     LLTROP,     I0+1,
     &               J0+1,      1,         CONV_CHK_SO2s )

         ! Update NSOFAR
         NSOFAR = NSOFAR + 2

!>>>
! Now include adjoint of F (dkh, 10/03/08) 

         ! Calculate NS (See DO_CONVECTION, NSTEP or NFCLDMX, NS
         CONVDT = GET_TS_CONV() * 60d0
         NSTEP  = CONVDT / 300
         NSTEP  = MAX( NSTEP, 1 )

         DO NS = 1, NSTEP
            ! QS_SO2 in NFCLDMX 
            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  NSOFAR,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,     JJPAR,     LLPAR,         I0+1,
     &                  J0+1,      1,         QC_SO2_CHK(:,:,:,NS) )

            ! Update NSOFAR
            NSOFAR = NSOFAR + 1

         ENDDO

         ! need this?  i don't think so...
         ! Update NSOFAR
         !NSOFAR = NSOFAR + 1
!<<<


      ENDIF ! LCONV
    
      ! Now checkpoint T_15_AVG and T_DAY so that we can use MEGAN emissions   (dkh, 01/22/10) 
      IF ( LMEGAN ) THEN 
     
         
         ! Only need to do this if it's a new day
         IF ( ITS_TIME_TO_CHK_T_15_AVG() ) THEN  

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO J = 1, JJPAR
            DO I = 1, IIPAR
 
               ! Get the values from megan_mod
               CHK_T_15_AVG(I,J,1) = (GET_T_15_AVG(I,J))

            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  1 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     1,             I0+1,
     &               J0+1,      1,         CHK_T_15_AVG )

            ! Set NSOFAR 
            NSOFAR = NSOFAR + 1
     
         ENDIF 

         ! Only need to do this if it's time for more A3
         IF ( ITS_TIME_FOR_A3_ADJ( ) ) THEN  

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, DAY_DIM
            DO J = 1, JJPAR
            DO I = 1, IIPAR
 
               ! Get the values from megan_mod
               CHK_T_DAY(I,J,L) = (GET_T_DAY(I,J,L))

            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO


            CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  1 + NSOFAR,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     DAY_DIM,       I0+1,
     &               J0+1,      1,         CHK_T_DAY )

            ! Set NSOFAR 
            NSOFAR = NSOFAR + 1
     
         ENDIF 
      ENDIF 

      ! Close file
      CLOSE( IU_RST )

      
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_CHK_DYN_FILE: wrote file' )

      ! Return to calling program
      END SUBROUTINE MAKE_CHK_DYN_FILE

!------------------------------------------------------------------------------

      SUBROUTINE READ_CHK_DYN_FILE( YYYYMMDD, HHMMSS ) 
!
!******************************************************************************
!  Subroutine READ_CHK_DYN_FILE reads values checkpointed at the dynamic time
!  step (dkh, 02/01/09) 
!
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day 
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!
!  Passed via CMN:
!  ============================================================================
!  (1 ) TAU    : TAU value (elapsed hours) at start of diagnostic interval
!
!  Notes
!  (1 ) Just like READ_CHK_DYN_FILE
!  (2 ) Add T_DAY and T_15_AVG (dkh, 01/23/10) 
!
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,         ONLY : OPEN_BPCH2_FOR_READ
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP, ALLOC_ERR
      USE FILE_MOD,          ONLY : IU_RST, IOERROR
      USE LOGICAL_MOD,       ONLY : LWETD, LCONV
      USE LOGICAL_MOD,       ONLY : LPRT
      USE LOGICAL_MOD,       ONLY : LMEGAN
      USE LOGICAL_ADJ_MOD,   ONLY : LAERO_THERM
      USE LOGICAL_ADJ_MOD,   ONLY : LDEL_CHKPT
      USE MEGAN_MOD,         ONLY : CHK_T_15_AVG
      USE MEGAN_MOD,         ONLY : CHK_T_DAY
      USE RESTART_MOD,       ONLY : CHECK_DIMENSIONS
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TIME_MOD,          ONLY : GET_TS_CONV
      USE TIME_MOD,          ONLY : ITS_TIME_TO_GET_T_15_AVG
      USE TIME_MOD,          ONLY : ITS_TIME_TO_GET_T_DAY
      USE TRACER_MOD,        ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,        ONLY : ITS_AN_AEROSOL_SIM
      USE UNIX_CMDS_MOD,     ONLY : REMOVE_CMD


#     include "CMN_SIZE"          ! Size parameters
#     include "comode.h"          ! ITLOOP, IGAS

      ! Arguments
      INTEGER, INTENT(IN)        :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER                    :: I, IOS, J, L, N, JLOOP, NN, NTL, AS
      INTEGER                    :: NCOUNT(NNPAR) 
      TYPE (XPLEX)                     :: TRACER(IIPAR,JJPAR,LLPAR)
      COMPLEX*16                     :: D_TRACER(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)                     :: SUMTC
      CHARACTER(LEN=255)         :: FILENAME
      CHARACTER(LEN=255)         :: UNZIP_FILE_CMD
      CHARACTER(LEN=255)         :: REMOVE_CHK_FILE_CMD
     

      ! For binary punch file, version 2.0
      INTEGER                    :: NI,     NJ,     NL
      INTEGER                    :: NS
      INTEGER                    :: NSTEP
      INTEGER                    :: CONVDT

      INTEGER                    :: IFIRST, JFIRST, LFIRST
      INTEGER                    :: NTRACER,   NSKIP
      INTEGER                    :: HALFPOLAR, CENTER180
      TYPE (XPLEX)                     :: LONRES,    LATRES
      TYPE (XPLEX)                     :: ZTAU0,     ZTAU1
      COMPLEX*16                     :: D_LONRES,    D_LATRES
      COMPLEX*16                     :: D_ZTAU0,     D_ZTAU1
      CHARACTER(LEN=20)          :: MODELNAME
      CHARACTER(LEN=40)          :: CATEGORY
      CHARACTER(LEN=40)          :: UNIT     
      CHARACTER(LEN=40)          :: RESERVED

      !=================================================================
      ! READ_CHECKPT_FILE begins here!
      !=================================================================

      INPUT_CHECKPT_FILE = 'gctm.chk.dyn.YYYYMMDD.hhmm'

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0
      D_TRACER(:,:,:) = 0e0
      !=================================================================
      ! Open checkpoint file and read top-of-file header
      !=================================================================
      
      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_CHECKPT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to name
      FILENAME = TRIM( ADJTMP_DIR ) // TRIM( FILENAME )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'C H E C K P T   F I L E   I N P U T'

      ! Remove obsolete options (dkh, 06/11/09) 
!      ! Unzip checkpt file
!      IF ( L_ZIP_CHECKPT ) THEN
!         UNZIP_FILE_CMD = TRIM( GUNZIP_CMD ) // ' ' //    
!     &                  TRIM( FILENAME ) // ZIP_SUFFIX    
!         CALL SYSTEM( TRIM( UNZIP_FILE_CMD ) )
!         WRITE( 6, 99 ) TRIM( UNZIP_FILE_CMD ) 
! 99   FORMAT( '     - READ_CHK_DYN_FILE: Executing: ',a )
!      ENDIF


      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_CHK_DYN_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
      
      !=================================================================
      ! Read checkpointed variables 
      !=================================================================
 
      ! Read in surface pressures
      READ( IU_RST, IOSTAT=IOS )
     &   MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
       LONRES%r=dble(D_LONRES)
       LATRES%r=dble(D_LATRES)
       LONRES%i=dimag(D_LONRES)
       LATRES%i=dimag(D_LATRES)
      ! IOS < 0 is end-of-file, so exit
      IF ( IOS < 0 ) GOTO 556

      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 )
     &   CALL IOERROR( IOS,IU_RST,'read_checkpt_file:10d' )

      READ( IU_RST, IOSTAT=IOS )
     &     CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &     NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &     NSKIP
       ZTAU0%r=dble(D_ZTAU0)
       ZTAU1%r=dble(D_ZTAU1)
       ZTAU0%i=dimag(D_ZTAU0)
       ZTAU1%i=dimag(D_ZTAU1)
      IF ( IOS /= 0 )
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:11d' )

      READ( IU_RST, IOSTAT=IOS )
     &     ( ( ( D_CHK_PSC(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
      CHK_PSC%r = dble(D_CHK_PSC)
      CHK_PSC%i = dimag(D_CHK_PSC)
      IF ( IOS /= 0 )
     &   CALL IOERROR( IOS,IU_RST,'read_checkpt_file:12d' )


      ! Read the values for WETDEP
      IF ( LWETD .and. 
     &   ( ITS_AN_AEROSOL_SIM() .or. ITS_A_FULLCHEM_SIM() ) ) THEN


         ! Read the values of WETD_CHK_H2O2s
         READ( IU_RST, IOSTAT=IOS )
     &   MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LONRES%r=dble(D_LONRES)
       LATRES%r=dble(D_LATRES)
       LONRES%i=dimag(D_LONRES)
       LATRES%i=dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 556

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:33d' )

         READ( IU_RST, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0%r=dble(D_ZTAU0)
       ZTAU1%r=dble(D_ZTAU1)
       ZTAU0%i=dimag(D_ZTAU0)
       ZTAU1%i=dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:34d' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( D_WETD_CHK_H2O2s(I,J,L), I=1,NI ), J=1,NJ ), 
     &              L=1,NL )
         WETD_CHK_H2O2s%r = dble(D_WETD_CHK_H2O2s)
         WETD_CHK_H2O2s%i = dimag(D_WETD_CHK_H2O2s)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:35d' )

         ! Read the values of WETD_CHK_SO2s
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LONRES%r=dble(D_LONRES)
       LATRES%r=dble(D_LATRES)
       LONRES%i=dimag(D_LONRES)
       LATRES%i=dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 556

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:36d' )

         READ( IU_RST, IOSTAT=IOS )
     &     CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0%r=dble(D_ZTAU0)
       ZTAU1%r=dble(D_ZTAU1)
       ZTAU0%i=dimag(D_ZTAU0)
       ZTAU1%i=dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:37d' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( D_WETD_CHK_SO2s(I,J,L), I=1,NI ), J=1,NJ ), 
     &               L=1,NL )
         WETD_CHK_SO2s%r=dble(D_WETD_CHK_SO2s)
         WETD_CHK_SO2s%i=dimag(D_WETD_CHK_SO2s)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:38d' )

         ! Read the values of WETD_CHK_SO4
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LONRES%r=dble(D_LONRES)
       LATRES%r=dble(D_LATRES)
       LONRES%i=dimag(D_LONRES)
       LATRES%i=dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 556

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:39d' )

         READ( IU_RST, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
         ZTAU0%r=dble(D_ZTAU0)
       ZTAU1%r=dble(D_ZTAU1)
       ZTAU0%i=dimag(D_ZTAU0)
       ZTAU1%i=dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:40d' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( D_WETD_CHK_SO4(I,J,L), I=1,NI ), J=1,NJ ), 
     &               L=1,NL )
         WETD_CHK_SO4%r=dble(D_WETD_CHK_SO4)
         WETD_CHK_SO4%i=dimag(D_WETD_CHK_SO4)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:41d' )

         ! Read the values of WETD_CHK_SO2
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LONRES%r=dble(D_LONRES)
       LATRES%r=dble(D_LATRES)
       LONRES%i=dimag(D_LONRES)
       LATRES%i=dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 556

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:42d' )

         READ( IU_RST, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
        ZTAU0%r=dble(D_ZTAU0)
       ZTAU1%r=dble(D_ZTAU1)
       ZTAU0%i=dimag(D_ZTAU0)
       ZTAU1%i=dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:43d' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( D_WETD_CHK_SO2(I,J,L), I=1,NI ), J=1,NJ ), 
     &               L=1,NL )
       WETD_CHK_SO2%r=dble(D_WETD_CHK_SO2)
       WETD_CHK_SO2%i=dimag(D_WETD_CHK_SO2)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:44d' )

      ENDIF ! LWETD

      ! Read the values for CONVECTION
      ! Replace NSRCX (dkh, 06/11/09) 
      !IF ( LCONV .AND. ( NSRCX == 3 .or. NSRCX == 10 ) ) THEN
      IF ( LCONV .AND. 
     &   ( ITS_AN_AEROSOL_SIM() .or. ITS_A_FULLCHEM_SIM() ) ) THEN


         ! Read the values of CONV_CHK_H2O2s
         READ( IU_RST, IOSTAT=IOS )
     &   MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LONRES%r=dble(D_LONRES)
       LATRES%r=dble(D_LATRES)
       LONRES%i=dimag(D_LONRES)
       LATRES%i=dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 556

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:51d' )

         READ( IU_RST, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
        ZTAU0%r=dble(D_ZTAU0)
       ZTAU1%r=dble(D_ZTAU1)
       ZTAU0%i=dimag(D_ZTAU0)
       ZTAU1%i=dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:52d' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( D_CONV_CHK_H2O2s(I,J,L), I=1,NI ), J=1,NJ ), 
     &              L=1,NL )
        CONV_CHK_H2O2s%r = dble(D_CONV_CHK_H2O2s)
        CONV_CHK_H2O2s%i = dimag(D_CONV_CHK_H2O2s)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:53d' )

         ! Read the values of CONV_CHK
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LONRES%r=dble(D_LONRES)
       LATRES%r=dble(D_LATRES)
       LONRES%i=dimag(D_LONRES)
       LATRES%i=dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 556

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:54d' )

         READ( IU_RST, IOSTAT=IOS )
     &     CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
        ZTAU0%r=dble(D_ZTAU0)
       ZTAU1%r=dble(D_ZTAU1)
       ZTAU0%i=dimag(D_ZTAU0)
       ZTAU1%i=dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:55d' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( D_CONV_CHK_SO2s(I,J,L), I=1,NI ), J=1,NJ ), 
     &               L=1,NL )
         CONV_CHK_SO2s%r = dble(D_CONV_CHK_SO2s)
         CONV_CHK_SO2s%i = dimag(D_CONV_CHK_SO2s)
         IF ( IOS /= 0 ) 
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:56d' )

!>>>
! Now include adjoint of F (dkh, 10/03/08) 

         ! Calculate NS (See DO_CONVECTION, NSTEP or NFCLDMX, NS
         CONVDT = GET_TS_CONV() * 60d0
         NSTEP  = CONVDT / 300
         NSTEP  = MAX( NSTEP, 1 )

         DO NS = 1, NSTEP

         ! Read the values of QC_SO2_CHK
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
         LONRES%r=dble(D_LONRES)
       LATRES%r=dble(D_LATRES)
       LONRES%i=dimag(D_LONRES)
       LATRES%i=dimag(D_LATRES)
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) GOTO 556

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 )
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:63d' )

         READ( IU_RST, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP
        ZTAU0%r=dble(D_ZTAU0)
       ZTAU1%r=dble(D_ZTAU1)
       ZTAU0%i=dimag(D_ZTAU0)
       ZTAU1%i=dimag(D_ZTAU1)
         IF ( IOS /= 0 )
     &      CALL IOERROR(IOS,IU_RST,'read_checkpt_file:64d' )

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( D_QC_SO2_CHK(I,J,L,NS), I=1,NI ), J=1,NJ ),
     &              L=1,NL )
        QC_SO2_CHK%r = dble(D_QC_SO2_CHK)
        QC_SO2_CHK%i = dimag(D_QC_SO2_CHK)
         IF ( IOS /= 0 )
     &      CALL IOERROR( IOS,IU_RST,'read_checkpt_file:65d' )


         ENDDO ! NS
!<<<


      ENDIF ! LCONV

      IF ( LMEGAN ) THEN
         
         ! adjoint equivalent of ITS_TIME_TO_CHK_T_15_AVG
         IF ( ITS_TIME_TO_GET_T_15_AVG() ) THEN 

            ! read in T_15_AVG
            READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
           LONRES%r=dble(D_LONRES)
       LATRES%r=dble(D_LATRES)
       LONRES%i=dimag(D_LONRES)
       LATRES%i=dimag(D_LATRES)
            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) GOTO 556

            ! IOS > 0 is a real I/O error -- print error message
            IF ( IOS > 0 ) 
     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:63' )

            READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP
          ZTAU0%r=dble(D_ZTAU0)
       ZTAU1%r=dble(D_ZTAU1)
       ZTAU0%i=dimag(D_ZTAU0)
       ZTAU1%i=dimag(D_ZTAU1)
            IF ( IOS /= 0 ) 
     &         CALL IOERROR(IOS,IU_RST,'read_checkpt_file:64' )

            READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), 
     &                 L=1,NL )
          TRACER%r = dble(D_TRACER)
          TRACER%i = dimag(D_TRACER)
            IF ( IOS /= 0 ) 
     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:59' )

            CHK_T_15_AVG(:,:,1) = TRACER(:,:,1)

         ENDIF 

         IF ( ITS_TIME_TO_GET_T_DAY( ) ) THEN 
    
            ! read in T_DAY
            READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
          LONRES%r=dble(D_LONRES)
       LATRES%r=dble(D_LATRES)
       LONRES%i=dimag(D_LONRES)
       LATRES%i=dimag(D_LATRES)
            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) GOTO 556

            ! IOS > 0 is a real I/O error -- print error message
            IF ( IOS > 0 ) 
     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:63' )

            READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP
          ZTAU0%r=dble(D_ZTAU0)
       ZTAU1%r=dble(D_ZTAU1)
       ZTAU0%i=dimag(D_ZTAU0)
       ZTAU1%i=dimag(D_ZTAU1)
            IF ( IOS /= 0 ) 
     &         CALL IOERROR(IOS,IU_RST,'read_checkpt_file:64' )

            READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), 
     &                 L=1,NL )
          TRACER%r = dble(D_TRACER)
          TRACER%i = dimag(D_TRACER)
            IF ( IOS /= 0 ) 
     &         CALL IOERROR( IOS,IU_RST,'read_checkpt_file:59' )

            CHK_T_DAY(:,:,:) = TRACER(:,:,:)

         ENDIF 

      ENDIF 
 556  CONTINUE 


      IF ( LDEL_CHKPT ) THEN 
         REMOVE_CHK_FILE_CMD  = TRIM ( REMOVE_CMD ) // ' ' //
     &                          TRIM ( FILENAME )

         CALL SYSTEM( TRIM( REMOVE_CHK_FILE_CMD ) )

         WRITE( 6, 102 ) TRIM( REMOVE_CHK_FILE_CMD )
 102     FORMAT( '     - READ_CHK_DYN_FILE: Executing: ',a )
      ENDIF 

      ! Remove obsolete (dkh, 06/11/09) 
!      ! Zip the .chk. file if it hasn't been deleted and zipping 
!      ! is requested
!      IF ( L_ZIP_CHECKPT .AND. (.NOT. L_DEL_CHECKPT) ) THEN
!         CALL BATCH_ZIP( YYYYMMDD, HHMMSS, 'chk', -1 )
!      ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_CHK_DYN_FILE: read file' )

      ! Return to calling program
      END SUBROUTINE READ_CHK_DYN_FILE

!-----------------------------------------------------------------------
      
      SUBROUTINE MAKE_ADJ_FILE( YYYYMMDD, HHMMSS )
!     
!******************************************************************************
!  Subroutine MAKE_ADJ_FILE creates a binary file of STT_ADJ
!  (dkh, 10/03/04)
!     
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Date
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create an adjoint file
!
!  Passed via CMN_ADJ 
!  ============================================================================
!  (1 ) ADJ_STT : Array of quantities to be checkpointed
!                               dim=(IIPAR,JJPAR,LLPAR,NADJ)
!
!  NOTES:
!  (1 ) Now write out adjoint of concentration scaling factors instead of 
!        adjoint of concentrations. This requires multiplying by STT.  This
!        routine is now called before chemistry and transport so that 
!        we can resale by the STT that was checkpointed after chemistry 
!        and transport in the forward run. STT is in [kg/box] at this point. 
!       Also, only write out LLADJKEEP levels and NNADJKEEP species
!        (dkh, 11/22/06)
!  (2 ) Update for GCv8 (dkh, 02/15/10) 
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : STT_ADJ
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJ_DIR
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,       ONLY : LPRT
      USE LOGICAL_ADJ_MOD,   ONLY : LTRAJ_SCALE
      USE TIME_MOD,          ONLY : EXPAND_DATE, GET_TAU
      USE TRACER_MOD,        ONLY : N_TRACERS 
      USE TRACER_MOD,        ONLY : STT

#     include "CMN_SIZE"       ! Size parameters
      
      ! Arguments
      INTEGER, INTENT(IN)     :: YYYYMMDD, HHMMSS
      
      ! Local Variables    
      INTEGER                 :: I,    I0, IOS, J,  J0, L, N
      INTEGER                 :: YYYY, MM, DD,  HH, SS
      TYPE (XPLEX)                  :: TRACER(IIPAR,JJPAR,LLPAR)
      CHARACTER(LEN=255)      :: FILENAME
      
      ! For binary punch file, version 2.0
      TYPE (XPLEX)                  :: LONRES, LATRES
      INTEGER, PARAMETER      :: HALFPOLAR = 1
      INTEGER, PARAMETER      :: CENTER180 = 1
      
      CHARACTER(LEN=40)       :: OUTPUT_ADJ_FILE
      CHARACTER(LEN=20)       :: MODELNAME
      CHARACTER(LEN=40)       :: CATEGORY
      CHARACTER(LEN=40)       :: UNIT
      CHARACTER(LEN=40)       :: RESERVED = ''
      CHARACTER(LEN=80)       :: TITLE
      
      ! Should make these user defined in input.gcadj
      !! Parameter
      INTEGER, PARAMETER      :: LLADJKEEP   = LLPAR
      !INTEGER, PARAMETER      :: NNADJKEEP   = N_TRACERS
       ! Now specify this input.gcadj
      !LOGICAL, PARAMETER      :: LTRAJ_SCALE = .TRUE. 
      
      !=================================================================
      ! MAKE_ADJ_FILE begins here!
      !=================================================================


      ! Hardwire output file for now
      OUTPUT_ADJ_FILE = 'gctm.adj.YYYYMMDD.hhmm'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM ADJ File: ' //
     &           'Instantaneous Adjoint Concentrations '
      CATEGORY = 'IJ-ADJ-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE
      
      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()
      
      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
      
      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================
      
      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_ADJ_FILE )
      
      ! Append the iteration number suffix to the file name
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )
      
      ! Add the ADJ_DIR prefix to the file name
      FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
      
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_ADJ_FILE: Writing ', a )
      
      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      !=================================================================
      ! Write each observed quantity to the observation file
      !=================================================================
      DO N = 1, N_TRACERS


         ! For saving out semilog sensitivities dJ/dSTT * STT = dJ/dln(STT)
         IF ( LTRAJ_SCALE ) THEN 
      
            UNIT     = 'J'

            !Temporarily store quantities in the TRACER array
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLADJKEEP
            DO J = 1, JJPAR
            DO I = 1, IIPAR
           
              ! Now multiply by concentrations so that we write 
              ! the adjoint of concentration scaling factors
              ! BUG FIX: it's better to use CHK_STT (dkh, 07/30/10) 
              !TRACER(I,J,L) = STT_ADJ (I,J,L,N) * STT(I,J,L,N)  
              TRACER(I,J,L) = STT_ADJ (I,J,L,N) * CHK_STT(I,J,L,N)  
            
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         
         ! For saving out sensitivities dJ/dSTT 
         ELSE 
      
            UNIT     = 'J/STT'

            !Temporarily store quantities in the TRACER array
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLADJKEEP
            DO J = 1, JJPAR
            DO I = 1, IIPAR
           
              ! Now multiply by concentrations so that we write 
              ! the adjoint of concentration scaling factors
              TRACER(I,J,L) = STT_ADJ (I,J,L,N) 
            
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         
         ENDIF 
 
         CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     LLADJKEEP, I0+1,
     &               J0+1,      1,         TRACER(:,:,1:LLADJKEEP)  )

      ENDDO
      
      ! Close file
      CLOSE( IU_RST )
      
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_ADJ_FILE: wrote file' )
      
      ! Return to calling program
      END SUBROUTINE MAKE_ADJ_FILE

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_EMS_ADJ_FILE( )
!
!******************************************************************************
!  Subroutine MAKE_EMS_ADJ_FILE creates a binary file of EMS_ADJ (dkh, 02/17/11) 
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC     : Current iteration number
!  (2 ) EMS_ADJ    : Array of adjoint gradients to be written
!
!  NOTES:
!  (1 ) Based on MAKE_GDT_FILE 
!
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,      ONLY : EMS_ADJ
      USE ADJ_ARRAYS_MOD,      ONLY : NNEMS, MMSCL 
      USE ADJ_ARRAYS_MOD,      ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,      ONLY : EXPAND_NAME
      USE ADJ_ARRAYS_MOD,      ONLY : COST_FUNC
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD,   ONLY : DIAGADJ_DIR
      USE ERROR_MOD,           ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,            ONLY : IU_RST,      IOERROR
      USE GRID_MOD,            ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_ADJ_MOD,     ONLY : LICS, LADJ_EMS 
      USE LOGICAL_MOD,         ONLY : LPRT
      USE TIME_MOD,            ONLY : EXPAND_DATE, GET_TAU
      USE TIME_MOD,            ONLY : GET_CT_EMIS

#     include "CMN_SIZE"   ! Size parameters
#     include "comode.h"   ! NEMIS(NCS)

      ! Local Variables
      INTEGER              :: I,    I0, IOS, J,  J0, L, M, N
      INTEGER              :: YYYY, MM, DD,  HH, SS
      CHARACTER(LEN=255)   :: FILENAME
     
      ! For binary punch file, version 2.0
      TYPE (XPLEX)               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      CHARACTER(LEN=20)    :: OUTPUT_FILE
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE

      !=================================================================
      ! MAKE_EMS_ADJ_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      OUTPUT_FILE = 'ems.adj.NN'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM ADJ File: ' //
     &           'Emissions adjoints '    
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_FILE )

      ! Append the iteration number suffix to the file name
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add the DIAGADJ_DIR prefix to the file name
      FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_EMS_ADJ_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      ! Set CATEGORY and UNIT
      CATEGORY = 'dJ_dEMS'
      UNIT     = 'J/kg'
 
      ! dkh debug
      !print*, ' CT_EMIS = ', GET_CT_EMIS()

      ! Convert units from J / (kg / box / timestep) to J / (kg / box)
      EMS_ADJ(:,:,:,:) = EMS_ADJ(:,:,:,:) / GET_CT_EMIS()

      !=================================================================
      ! Write each observed quantity to the observation file
      !=================================================================
      DO N = 1, NNEMS

         CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     MMSCL,     I0+1,
     &               J0+1,      1,     (EMS_ADJ(:,:,:,N)) )

      ENDDO

      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_EMS_ADJ_FILE: wrote file' )

      ! Return to calling program
      END SUBROUTINE MAKE_EMS_ADJ_FILE

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_PROD_GDT_FILE( )
!           
!******************************************************************************
!  Subroutine MAKE_PROD_GDT_FILE (GDT=SF_ADJ) creates a binary file of 
!  PROD_SF_ADJ (hml, 07/26/11, adj32_025)
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC     : Current iteration number
!  (2 ) PROD_ADJ    : Array of adjoint gradients to be written
!           
!  NOTES:   
!  (1 ) Based on MAKE_EMS_ADJ_FILE & MAKE_ADJ_FILE 
!
!******************************************************************************
!        
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,      ONLY : PROD_SF_ADJ
      USE ADJ_ARRAYS_MOD,      ONLY : NSTPL, MMSCL
      USE ADJ_ARRAYS_MOD,      ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,      ONLY : EXPAND_NAME
      USE ADJ_ARRAYS_MOD,      ONLY : COST_FUNC
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD,   ONLY : DIAGADJ_DIR
      USE ERROR_MOD,           ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,            ONLY : IU_RST,      IOERROR
      USE GRID_MOD,            ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,         ONLY : LPRT
      USE TIME_MOD,            ONLY : EXPAND_DATE, GET_TAU

#     include "CMN_SIZE"   ! Size parameters
#     include "comode.h"   ! NEMIS(NCS)

      ! Local Variables
      INTEGER              :: I,    I0, IOS, J,  J0, L, M, N
      INTEGER              :: YYYY, MM, DD,  HH, SS
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      TYPE (XPLEX)               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      CHARACTER(LEN=40)    :: OUTPUT_FILE
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE


      !=================================================================
      ! MAKE_PROD_ADJ_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      OUTPUT_FILE = 'prod.sf.adj.NN'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM ADJ File: ' //
     &           'Stratospheric Production Scaling Factor Adjoints '
      UNIT     = 'J'
      !CATEGORY = 'dJ_dPRSF'
      CATEGORY = 'IJ-GDP-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_FILE )

      ! Append the iteration number suffix to the file name
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add the DIAGADJ_DIR prefix to the file name
      FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_PROD_GDT_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      !=================================================================
      ! Write each observed quantity to the observation file
      !=================================================================
      DO N = 1, NSTPL

         CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     MMSCL,     I0+1,
     &               J0+1,      1,   (PROD_SF_ADJ(:,:,1,N)))

      ENDDO

      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG(
     &                 '### MAKE_PROD_GDT_FILE: wrote file' )

      ! Return to calling program
      END SUBROUTINE MAKE_PROD_GDT_FILE

!------------------------------------------------------------------------------
!
      SUBROUTINE MAKE_LOSS_GDT_FILE( )
!           
!******************************************************************************
!  Subroutine MAKE_LOSS_GDT_FILE (GDT = SF_ADJ) creates a binary file of 
!  LOSS_SF_ADJ: stratospheric loss rate scaling factor adjoint
!  (hml, 07/26/11, adj32_025)
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC     : Current iteration number
!  (2 ) LOSS_ADJ    : Array of adjoint gradients to be written
!           
!           
!  NOTES:   
!  (1 ) Based on MAKE_EMS_ADJ_FILE & MAKE_ADJ_FILE 
!
!******************************************************************************
!        
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,      ONLY : LOSS_SF_ADJ
      USE ADJ_ARRAYS_MOD,      ONLY : NSTPL, MMSCL
      USE ADJ_ARRAYS_MOD,      ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,      ONLY : EXPAND_NAME
      USE ADJ_ARRAYS_MOD,      ONLY : COST_FUNC
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD,   ONLY : DIAGADJ_DIR
      USE ERROR_MOD,           ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,            ONLY : IU_RST,      IOERROR
      USE GRID_MOD,            ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,         ONLY : LPRT
      USE TIME_MOD,            ONLY : EXPAND_DATE, GET_TAU

#     include "CMN_SIZE"   ! Size parameters
#     include "comode.h"   ! NEMIS(NCS)

      ! Local Variables
      INTEGER              :: I,    I0, IOS, J,  J0, L, M, N
      INTEGER              :: YYYY, MM, DD,  HH, SS
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      TYPE (XPLEX)               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      CHARACTER(LEN=40)    :: OUTPUT_FILE
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE

      !=================================================================
      ! MAKE_LOSS_SF_ADJ_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      OUTPUT_FILE = 'loss.sf.adj.NN'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM ADJ File: ' //
     &           'Stratospheric Loss Scaling Factor Adjoints '
      !CATEGORY = 'dJ_dLSSF'
      CATEGORY = 'IJ-GDL-$'
      UNIT     = 'J'

      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_FILE )

      ! Append the iteration number suffix to the file name
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add the DIAGADJ_DIR prefix to the file name
      FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_LOSS_GDT_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      !=================================================================
      ! Write each observed quantity to the observation file
      !=================================================================
      DO N = 1, NSTPL

         CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     MMSCL,     I0+1,
     &               J0+1,      1,    (LOSS_SF_ADJ(:,:,1,N)))

      ENDDO

      ! Close file
      CLOSE( IU_RST )

      !### Debug 
      IF ( LPRT ) CALL DEBUG_MSG(
     &                 '### MAKE_LOSS_GDT_FILE: wrote file' )
      
      ! Return to calling program
      END SUBROUTINE MAKE_LOSS_GDT_FILE

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_SO2ac_FILE( ESO2_ac, MONTH )
!
!******************************************************************************
!  Subroutine MAKE_SO2ac_FILE creates GEOS-CHEM checkpt files of SO2 
!  emissions from aircraft. 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ESO2_ac  :  Current monthly aircraft SO2 emissions [kg SO2/box/s]
!  (2 ) MONTH    :  Current month 
!
!  NOTES:
!
!******************************************************************************
!     
      ! References to F90 modules
      USE BPCH2_MOD
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP, ALLOC_ERR
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR
      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,       ONLY : LPRT 
      USE TIME_MOD,          ONLY : EXPAND_DATE, GET_TAU
      USE TIME_MOD,          ONLY : GET_NYMD, GET_NHMS 

#     include "CMN_SIZE"          ! Size parameters

      ! Arguments
      CHARACTER(LEN=3), INTENT(IN) :: MONTH
      TYPE (XPLEX),           INTENT(IN) :: ESO2_ac(IIPAR,JJPAR,LLPAR)

      ! Local Variables      
      INTEGER              :: I,    I0, IOS, J,  J0, L
      INTEGER              :: YYYY, MM, DD,  HH, SS, ZIP_HH, AS
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      TYPE (XPLEX)               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT     
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE 

      !=================================================================
      ! MAKE_SO2ac_FILE begins here!
      !=================================================================

      ! NEW: rename them *.chk.con.*  (dkh, 10/10/08) 
      OUTPUT_CHECKPT_FILE = 'gctm.chk.SO2ac.YYYY.' // TRIM( MONTH ) 

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM SO2 ac     Checkpoint File '
      CATEGORY = 'IJ-CHK-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the checkpoint file for output -- binary punch format
      !=================================================================

      ! Copy the output checkpoint file name into a local variable
      FILENAME = TRIM( OUTPUT_CHECKPT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )

      ! Add ADJ_DIR prefix to filename
      FILENAME = TRIM( ADJTMP_DIR ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_SO2ac_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      !=================================================================
      ! Write each checkpointed quantity to the checkpoint file
      !=================================================================
    
      CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  1,
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     LLPAR,     I0+1,
     &            J0+1,      1,     ESO2_ac )


      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_SO2ac_FILE: wrote file' )

      ! Return to calling program
      END SUBROUTINE MAKE_SO2ac_FILE

!------------------------------------------------------------------------------
      SUBROUTINE READ_SO2ac_FILE( ESO2_ac, MONTH )
!
!******************************************************************************
!  Subroutine READ_SO2ac_FILE initializes GEOS-CHEM tracer concentrations 
!  from a checkpoint file (binary punch file format) from before convection
!  (dkh, 8/30/04, mak, 8/2/07)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) MONTH    :  Current month
!
!  Arguments as output:
!  ============================================================================
!  (1 ) ESO2_ac  :  Current monthly aircraft SO2 emissions [kg SO2/box/s]
!
!  Notes
!
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,         ONLY : OPEN_BPCH2_FOR_READ
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR 
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP, ALLOC_ERR
      USE FILE_MOD,          ONLY : IU_RST, IOERROR
      USE LOGICAL_MOD,       ONLY : LPRT
      USE LOGICAL_ADJ_MOD,   ONLY : LDEL_CHKPT
      USE RESTART_MOD,       ONLY : CHECK_DIMENSIONS
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TIME_MOD,          ONLY : GET_NYMD
      USE TIME_MOD,          ONLY : GET_NHMS
      USE UNIX_CMDS_MOD,     ONLY : REMOVE_CMD


#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      CHARACTER(LEN=3), INTENT(IN) :: MONTH
      TYPE (XPLEX),          INTENT(OUT) :: ESO2_ac(IIPAR,JJPAR,LLPAR)

      ! Local Variables
      INTEGER             :: I, IOS, J, L, N, JLOOP, NN, NTL, AS
      TYPE (XPLEX)              :: TRACER(IIPAR,JJPAR,LLPAR)
      COMPLEX*16              :: D_TRACER(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME
      CHARACTER(LEN=255)  :: REMOVE_CHK_FILE_CMD
     

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      TYPE (XPLEX)              :: LONRES,    LATRES
      TYPE (XPLEX)              :: ZTAU0,     ZTAU1
      COMPLEX*16              :: D_LONRES,    D_LATRES
      COMPLEX*16              :: D_ZTAU0,     D_ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     
      CHARACTER(LEN=40)   :: RESERVED

      !=================================================================
      ! READ_SO2ac_FILE begins here!
      !=================================================================

      INPUT_CHECKPT_FILE = 'gctm.chk.SO2ac.YYYY.' // TRIM( MONTH ) 

      ! Initialize some variables
      TRACER(:,:,:) = 0e0

      !=================================================================
      ! Open checkpoint file and read top-of-file header
      !=================================================================
      
      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_CHECKPT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )

      ! Add ADJ_DIR prefix to name
      FILENAME = TRIM( ADJTMP_DIR ) // TRIM( FILENAME )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'C H E C K P T   F I L E   I N P U T'

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_SO2ac: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )
      
      !=================================================================
      ! Read checkpointed variables 
      !=================================================================
 
      READ( IU_RST, IOSTAT=IOS )
     &  MODELNAME, D_LONRES, D_LATRES, HALFPOLAR, CENTER180
       LATRES%r=dble(D_LATRES)
       LONRES%r=dble(D_LONRES)
       LATRES%i=dimag(D_LATRES)
       LONRES%i=dimag(D_LONRES)
      ! IOS < 0 is end-of-file, so print error message
      IF ( IOS < 0 ) 
     &   CALL IOERROR( IOS,IU_RST,'read_so2ac_file:6' )

      ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 ) 
     &   CALL IOERROR( IOS,IU_RST,'read_so2ac_file:7' )

      READ( IU_RST, IOSTAT=IOS )
     &     CATEGORY, NTRACER,  UNIT, D_ZTAU0,  D_ZTAU1,  RESERVED,
     &     NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &     NSKIP
        ZTAU0%r=dble(D_ZTAU0)
       ZTAU1%r=dble(D_ZTAU1)
       ZTAU0%i=dimag(D_ZTAU0)
       ZTAU1%i=dimag(D_ZTAU1)
         IF ( IOS /= 0 ) 
     &   CALL IOERROR(IOS,IU_RST,'read_so2ac_file:8' )

      READ( IU_RST, IOSTAT=IOS )
     &     ( ( ( D_TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
       TRACER%r=dble(D_TRACER)
     
       TRACER%i=dimag(D_TRACER)
    
      IF ( IOS /= 0 ) 
     &   CALL IOERROR( IOS,IU_RST,'read_so2ac_file:9' )

      !==============================================================
      ! Assign data from the TRACER array to the STT array.
      !==============================================================

      ! Only process checkpoint data (i.e. mixing ratio)
      IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN

         ! Make sure array dimensions are of global size
         ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
         !print*, 'before check_dimensions ni, nj, nl are', ni, nj, nl
         CALL CHECK_DIMENSIONS( NI, NJ, NL )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ESO2_ac(I,J,L) = TRACER(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF !category is checkpoint

      ! Close file
      CLOSE( IU_RST )      


 555  CONTINUE

      ! Remove files if L_CHK_DEL = TRUE 
      IF ( LDEL_CHKPT ) THEN 
         REMOVE_CHK_FILE_CMD  = TRIM ( REMOVE_CMD ) // ' ' //
     &                          TRIM ( FILENAME )

        CALL SYSTEM( TRIM( REMOVE_CHK_FILE_CMD ) )

        WRITE( 6, 102 ) TRIM( REMOVE_CHK_FILE_CMD )
 102    FORMAT( '     - READ_SO2ac_FILE: Executing: ',a )
      ENDIF 


      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_SO2ac_FILE: read file' )

      ! Return to calling program
      END SUBROUTINE READ_SO2ac_FILE

!-----------------------------------------------------------------------
      SUBROUTINE EXPAND_NAME( FILENAME, N_ITRN )
!
!******************************************************************************
!  Subroutine EXPAND_DATE replaces "NN" token within
!  a filename string with the actual values. (bmy, 6/27/02, 12/2/03)
!  (dkh, 9/22/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Filename with tokens to replace
!  (2 ) N_ITRN   (INTEGER  ) : Current iteration number
!
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Modified filename
!
!  NOTES:
!  (1 ) Based on EXPAND_DATE
!
!******************************************************************************
!
      ! References to F90 modules
      USE CHARPAK_MOD, ONLY : STRREPL
      USE ERROR_MOD,   ONLY : ERROR_STOP

#     include "define.h"

      ! Arguments
      CHARACTER(LEN=*), INTENT(INOUT) :: FILENAME
      INTEGER,          INTENT(IN)    :: N_ITRN

      ! Local variables
      CHARACTER(LEN=2)                :: NN_STR

      !=================================================================
      ! EXPAND_NAME begins here!
      !=================================================================

#if   defined( LINUX_PGI )

      ! Use ENCODE statement for PGI/Linux (bmy, 9/29/03)
      ENCODE( 2, '(i2.2)', NN_STR   ) N_ITRN

#else

      ! For other platforms, use an F90 internal write (bmy, 9/29/03)
      WRITE( NN_STR,   '(i2.2)' ) N_ITRN

#endif

      ! Replace NN token w/ actual value
      CALL STRREPL( FILENAME, 'NN',   NN_STR   )


      ! Return to calling program
      END SUBROUTINE EXPAND_NAME

!-----------------------------------------------------------------------
      SUBROUTINE INIT_CHECKPT
!
!*****************************************************************************
!  Subroutine INIT_CHECKPT initializes all module arrays (dkh, 9/10/04)
!
!  NOTES:
!  (1 ) Add CHK_PSC. (dkh, 03/16/05)
!  (2 ) Add ORIG_STT. (dkh, 06/14/05)
!  (3 ) Add PART_CASE. (dkh, 07/22/05)
!  (4 ) Add CHK_STT_BEFCHEM. (dkh, 08/08/05)
!  (5 ) Add SO2_CHK, H2O2_CHK. (dkh, 10/23/05) 
!  (6 ) Add WETD_CHK_H2O2s_CHEMT, WETD_CHK_H2O2s_DYNT, etc. (dkh, 10/23/05) 
!  (7 ) Add WETD_CHK_SO2_CHEMT, WETD_CHK_SO2_DYNT. (dkh, 10/31/05)  
!  (8 ) Add CONV_CHK... (dkh, 11/22/05)  
!  (9 ) Add SOILNOX_CHK (dkh, 02/06/07) 
!  (10) Change to checkpointing WETD and CONV stuff at every dynamic ts (dkh, 02/02/09) 
!  (11) Update to v8, (adj_group, 6/09/09)
!  (12) Add CHK_T_DAY and CHK_T_15_DAY for MEGAN emissions (dkh, 01/23/10) 
!******************************************************************************
!      
      ! F90 modules
      USE ERROR_MOD,        ONLY : ALLOC_ERR
      USE TIME_MOD,         ONLY : GET_TS_CONV
      !USE ADJ_ARRAYS_MOD,   ONLY : NOBS
      USE TRACER_MOD,       ONLY : N_TRACERS
      USE TRACER_MOD,       ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,       ONLY : ITS_AN_AEROSOL_SIM
      USE LOGICAL_MOD,      ONLY : LSOILNOX
      USE LOGICAL_MOD,      ONLY : LSULF
      USE LOGICAL_MOD,      ONLY : LCHEM
      USE LOGICAL_MOD,      ONLY : LWETD, LCONV
      USE LOGICAL_MOD,      ONLY : LMEGAN
      USE LOGICAL_MOD,      ONLY : LLIGHTNOX
      USE LOGICAL_MOD,      ONLY : LWETD
      USE LOGICAL_ADJ_MOD,  ONLY : LAERO_THERM
      USE MEGAN_MOD,        ONLY : DAY_DIM


#     include "CMN_SIZE"         ! IIPAR, JJPAR, LLPAR
#     include "comode.h"         ! ITLOOP

      ! Local variables
      INTEGER                   :: AS      
      INTEGER                   :: NSTEP
      INTEGER                   :: CONVDT


      !=================================================================      
      ! INIT_CHECKPT begins here
      !=================================================================      
     
      IF ( LSULF .and. LAERO_THERM ) THEN 
         ALLOCATE( RP_IN( IIPAR, JJPAR, LLPAR, NRPIN ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'RP_IN' )

         ALLOCATE( RP_OUT( IIPAR, JJPAR, LLPAR, NRPOUT ) , STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'RP_OUT' )

         ALLOCATE( nitr_max( IIPAR, JJPAR, LLPAR ) , STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'nitr_max' )

         ALLOCATE( gamaan_fwd( IIPAR, JJPAR, LLPAR, NNNMAX ) , STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'gamaan_fwd' )

         ALLOCATE( gamold_fwd( IIPAR, JJPAR, LLPAR, NNNMAX ) , STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'gamold_fwd' )

         ALLOCATE( wh2o_fwd( IIPAR, JJPAR, LLPAR, NNNMAX ) , STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'wh2o_fwd' )
   
         ALLOCATE( ynh4_fwd( IIPAR, JJPAR, LLPAR, NNNMAX ) , STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ynh4_fwd' )
   
         ALLOCATE( eror_fwd( IIPAR, JJPAR, LLPAR, NNNMAX ) , STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'error_fwd' )
   
         ALLOCATE( exit_fwd( IIPAR, JJPAR, LLPAR, NNNMAX ) , STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'exit_fwd' )
   
         ALLOCATE( gamana_fwd( IIPAR, JJPAR, LLPAR, NNNMAX ) , STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'gaman_fwd' )
   
         ALLOCATE( gamas1_fwd( IIPAR, JJPAR, LLPAR, NNNMAX ) , STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'gamas1_fwd' )
   
         ALLOCATE( gamas2_fwd( IIPAR, JJPAR, LLPAR, NNNMAX ) , STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'gamas2_fwd' )

      ENDIF 
 
      IF ( LSULF .and. LCHEM ) THEN 
         ALLOCATE( SO2_CHK( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2_CHK' )
         SO2_CHK = 0.d0
         ALLOCATE( D_SO2_CHK( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'D_SO2_CHK' )
         D_SO2_CHK = 0.d0 

         ALLOCATE( H2O2_CHK( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'H2O2_CHK' )
         H2O2_CHK = 0.d0
         ALLOCATE( D_H2O2_CHK( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'D_H2O2_CHK' )
         D_H2O2_CHK = 0.d0 

      ENDIF 

      ALLOCATE( CHK_STT( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CHK_STT' )

      ! mak
      ALLOCATE( CHK_STT_CON( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CHK_STT_CON' )

      ! OBS_STT now in adj_arrays_mod.f instead of checkpt_mod.f (mak, 6/14/09)
      !ALLOCATE( OBS_STT( IIPAR, JJPAR, LLPAR, N_TRACERS ) , STAT=AS )
      !IF ( AS /= 0 ) CALL ALLOC_ERR( 'OBS' )

      ALLOCATE( CHK_PSC( IIPAR, JJPAR, 2 ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CHK_PSC' )

      ALLOCATE( D_CHK_PSC( IIPAR, JJPAR, 2 ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'D_CHK_PSC' )

      ALLOCATE( ORIG_STT( IIPAR, JJPAR, LLPAR, N_TRACERS ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ORIG_STT' )

      IF ( ITS_A_FULLCHEM_SIM() .and. 
     &     ( LCHEM .or. LWETD )        ) THEN 
         ALLOCATE( PART_CASE( ITLOOP ) , STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PART_CASE' )

         ALLOCATE( CHK_STT_BEFCHEM( IIPAR, JJPAR, LLPAR, N_TRACERS ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CHK_STT_BEFCHEM' )

         ALLOCATE( CHK_HSAVE( IIPAR, JJPAR, LLTROP ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CHK_HSAVE' )
         CHK_HSAVE = 0.d0
         ALLOCATE( D_CHK_HSAVE( IIPAR, JJPAR, LLTROP ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'D_CHK_HSAVE' )
         D_CHK_HSAVE = 0.d0 
         
      ENDIF 

      IF ( LWETD .and. 
     &   ( ITS_AN_AEROSOL_SIM() .or. ITS_A_FULLCHEM_SIM() ) ) THEN 
         ALLOCATE( WETD_CHK_H2O2s( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'H2O2s' )
         ALLOCATE( D_WETD_CHK_H2O2s( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'H2O2s' )
         D_WETD_CHK_H2O2s = 0.d0
         WETD_CHK_H2O2s = 0.d0 

         ALLOCATE( WETD_CHK_SO2s( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2s' )
         WETD_CHK_SO2s = 0.d0
         ALLOCATE( D_WETD_CHK_SO2s( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2s' )
         D_WETD_CHK_SO2s = 0.d0 

         ALLOCATE( WETD_CHK_SO4( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO4' )
         WETD_CHK_SO4 = 0.d0 
         ALLOCATE( D_WETD_CHK_SO4( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO4' )
         D_WETD_CHK_SO4 = 0.d0

         ALLOCATE( WETD_CHK_SO2( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2' )
         WETD_CHK_SO2 = 0.d0
         ALLOCATE( D_WETD_CHK_SO2( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2' )
         D_WETD_CHK_SO2 = 0.d0 

      ENDIF  ! LWETD 

      IF ( LCONV .AND. 
     &   ( ITS_AN_AEROSOL_SIM() .or. ITS_A_FULLCHEM_SIM() ) ) THEN 

         ALLOCATE( CONV_CHK_H2O2s( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CONV_H2O2s' )
         CONV_CHK_H2O2s = 0.d0 
         ALLOCATE( D_CONV_CHK_H2O2s( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CONV_H2O2s' )
         D_CONV_CHK_H2O2s = 0.d0


         ALLOCATE( CONV_CHK_SO2s( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CONV_SO2s' )
         CONV_CHK_SO2s = 0.d0 

         ALLOCATE( D_CONV_CHK_SO2s( IIPAR, JJPAR, LLPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CONV_SO2s' )
         D_CONV_CHK_SO2s = 0.d0
      !>>>                
      ! Now include adjoint of F (dkh, 10/03/08) 

         ! Calculate NS (See DO_CONVECTION, NSTEP or NFCLDMX, NS
         CONVDT = GET_TS_CONV() * 60d0
         NSTEP  = CONVDT / 300
         NSTEP  = MAX( NSTEP, 1 )

         ALLOCATE( QC_SO2_CHK( IIPAR, JJPAR, LLPAR, NSTEP),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'QC_SO2_CHK' )
         QC_SO2_CHK = 0.d0
         ALLOCATE( D_QC_SO2_CHK( IIPAR, JJPAR, LLPAR, NSTEP),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'D_QC_SO2_CHK' )
         D_QC_SO2_CHK = 0.d0

      !<<< 

      ENDIF  ! LCONV

      IF ( LSOILNOX ) THEN

         ALLOCATE( SOILNOX_CHK( IIPAR, JJPAR ),
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'SOILNOX_CHK' )
         SOILNOX_CHK = 0.d0 

      ENDIF

      ! Adding this didn't really help (dkh, 06/11/09) 
      !IF ( LADJ_TRAN ) THEN
      !
      !   ALLOCATE( CHK_STT_TD( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      !   IF ( AS /= 0 ) CALL ALLOC_ERR( 'CHK_STT_TD' )
      !   CHK_STT_TD = 0.d0
      !
      !   ALLOCATE( CHK_STT_TC( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      !   IF ( AS /= 0 ) CALL ALLOC_ERR( 'CHK_STT_TC' )
      !   CHK_STT_TC = 0.d0
      !
      !ENDIF

      ! adj_group: add for checkpointing lightning NOx emissions
      IF ( LLIGHTNOX ) THEN 
         ALLOCATE( SLBASE_CHK( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'SLBASE_CHK' )
         SLBASE_CHK = 0d0
      ENDIF 

      END SUBROUTINE INIT_CHECKPT

!-----------------------------------------------------------------------

      SUBROUTINE CLEANUP_CHECKPT
!
!*****************************************************************************
!  Subroutine CLEANUP_CHECKPT deallocates all module arrays (dkh, 9/10/04)
!
!  NOTES:
!  (1 ) Add CHK_PSC. (dkh, 03/16/05) 
!  (2 ) Add ORIG_STT. (dkh, 06/14/05)
!  (3 ) Add PART_CASE. (dkh, 07/22/05)
!  (4 ) Add CHK_STT_BEFCHEM. (dkh, 08/08/05) 
!  (5 ) Add SO2_CHK, H2O2_CHK. (dkh, 10/23/05) 
!  (6 ) Add WETD_CHK_H2O2s_CHEMT, WETD_CHK_H2O2s_DYNT, etc. (dkh, 10/23/05) 
!  (7 ) Add WETD_CHK_SO2_CHEMT, WETD_CHK_SO2_DYNT. (dkh, 10/31/05)  
!  (8 ) Add CONV_CHK_xxx etc.  (dkh, 11/22/05)  
!  (9 ) Add SOILNOX_CHK. (dkh, 02/06/07) 
!  (10) Change to checkpointing WETD and CONV stuff at every dynamic ts (dkh, 02/02/09) 
!******************************************************************************
!
      IF ( ALLOCATED( RP_IN ) )      DEALLOCATE( RP_IN )
      IF ( ALLOCATED( RP_OUT) )      DEALLOCATE( RP_OUT )
      !OBS_STT now in adj_arrays_mod.f instead of checkpt_mod.f (mak, 6/14/09)
      !IF ( ALLOCATED( OBS_STT ) )    DEALLOCATE( OBS_STT )
      IF ( ALLOCATED( CHK_STT ) )    DEALLOCATE( CHK_STT )
      IF ( ALLOCATED( CHK_PSC ) )    DEALLOCATE( CHK_PSC )
      IF ( ALLOCATED( D_CHK_PSC ) )    DEALLOCATE( D_CHK_PSC )
      IF ( ALLOCATED( nitr_max ) )   DEALLOCATE( nitr_max )
      IF ( ALLOCATED( gamaan_fwd ) ) DEALLOCATE( gamaan_fwd )
      IF ( ALLOCATED( gamold_fwd ) ) DEALLOCATE( gamold_fwd )
      IF ( ALLOCATED( wh2o_fwd ) )   DEALLOCATE( wh2o_fwd )
      IF ( ALLOCATED( ynh4_fwd ) )   DEALLOCATE( ynh4_fwd )
      IF ( ALLOCATED( eror_fwd ) )   DEALLOCATE( eror_fwd )
      IF ( ALLOCATED( exit_fwd ) )   DEALLOCATE( exit_fwd )
      IF ( ALLOCATED( gamana_fwd ) ) DEALLOCATE( gamana_fwd )
      IF ( ALLOCATED( gamas1_fwd ) ) DEALLOCATE( gamas1_fwd )
      IF ( ALLOCATED( gamas2_fwd ) ) DEALLOCATE( gamas2_fwd )
      IF ( ALLOCATED( ORIG_STT ) )   DEALLOCATE( ORIG_STT )
      IF ( ALLOCATED( CHK_STT_BEFCHEM ) )  DEALLOCATE( CHK_STT_BEFCHEM )
      IF ( ALLOCATED( PART_CASE ) )  DEALLOCATE( PART_CASE )
      IF ( ALLOCATED( CHK_HSAVE ) )  DEALLOCATE( CHK_HSAVE )
      IF ( ALLOCATED( D_CHK_HSAVE ) )  DEALLOCATE( D_CHK_HSAVE )
      IF ( ALLOCATED( SO2_CHK ) )  DEALLOCATE( SO2_CHK )
      IF ( ALLOCATED( D_SO2_CHK ) )  DEALLOCATE( D_SO2_CHK )
      IF ( ALLOCATED( H2O2_CHK ) )  DEALLOCATE( H2O2_CHK )
      IF ( ALLOCATED( D_H2O2_CHK ) )  DEALLOCATE( D_H2O2_CHK )
      IF ( ALLOCATED( WETD_CHK_H2O2s ) )  
     &    DEALLOCATE( WETD_CHK_H2O2s )
      IF ( ALLOCATED( D_WETD_CHK_H2O2s ) )
     &    DEALLOCATE( D_WETD_CHK_H2O2s )
      IF ( ALLOCATED( WETD_CHK_SO2s ) )  
     &    DEALLOCATE( WETD_CHK_SO2s )
      IF ( ALLOCATED( D_WETD_CHK_SO2s ) )
     &    DEALLOCATE( D_WETD_CHK_SO2s )
      IF ( ALLOCATED( WETD_CHK_SO4 ) )  
     &    DEALLOCATE( WETD_CHK_SO4 )
      IF ( ALLOCATED( D_WETD_CHK_SO4 ) )
     &    DEALLOCATE( D_WETD_CHK_SO4 )
      IF ( ALLOCATED( WETD_CHK_SO2 ) )  
     &    DEALLOCATE( WETD_CHK_SO2 )
      IF ( ALLOCATED( D_WETD_CHK_SO2 ) )
     &    DEALLOCATE( D_WETD_CHK_SO2 )
      IF ( ALLOCATED( CONV_CHK_H2O2s ) )  
     &    DEALLOCATE( CONV_CHK_H2O2s )
      IF ( ALLOCATED( D_CONV_CHK_H2O2s ) )
     &    DEALLOCATE( D_CONV_CHK_H2O2s )
      IF ( ALLOCATED( CONV_CHK_SO2s ) )  
     &    DEALLOCATE( CONV_CHK_SO2s )
      IF ( ALLOCATED( D_CONV_CHK_SO2s ) )
     &    DEALLOCATE( D_CONV_CHK_SO2s )
      IF ( ALLOCATED( SOILNOX_CHK ) )  
     &    DEALLOCATE( SOILNOX_CHK )

      IF ( ALLOCATED( CHK_STT_CON ) )DEALLOCATE( CHK_STT_CON )

      !IF ( ALLOCATED( CHK_STT_TD ) )DEALLOCATE( CHK_STT_TD )
      !IF ( ALLOCATED( CHK_STT_TC ) )DEALLOCATE( CHK_STT_TC )

!>>>            
! Now include adjoint of F (dkh, 10/03/08) 
      IF ( ALLOCATED( QC_SO2_CHK) )
     &    DEALLOCATE( QC_SO2_CHK )
      IF ( ALLOCATED( D_QC_SO2_CHK) )
     &    DEALLOCATE( D_QC_SO2_CHK )
!<<<


      IF ( ALLOCATED( SLBASE_CHK   ) )  DEALLOCATE( SLBASE_CHK   )


      ! Return to calling program
      END SUBROUTINE CLEANUP_CHECKPT

!------------------------------------------------------------------------------
      END MODULE CHECKPT_MOD

