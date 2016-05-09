! $Id: upbdflx_adj_mod.f,v 1.3 2012/03/01 22:00:26 daven Exp $
      MODULE UPBDFLX_ADJ_MOD
!
!******************************************************************************
!  Module UPBDFLX_MOD contains subroutines which impose stratospheric boundary
!  conditions on O3 and NOy (qli, bdf, mje, bmy, 6/28/01, 12/1/04)
!
!  Module Variables:
!  ===========================================================================
!  (1 ) IORD (INTEGER) : TPCORE E/W      transport option flag 
!  (2 ) JORD (INTEGER) : TPCORE N/S      transport option flag
!  (3 ) KORD (INTEGER) : TPCORE vertical transport option flag
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_UPBDFLX     : Driver for stratospheric flux boundary conditions
!  (2 ) UPBDFLX_O3     : Computes flux of O3 from stratosphere, using Synoz
!  (3 ) UPBDFLX_NOY    : Computes flux of NOy from stratosphere
!  (4 ) INIT_UPBDFLX   : Gets IORD, JORD, KORD values from "input_mod.f"
!
!  GEOS-CHEM modules referenced by upbdflx_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f    : Module containing routines for binary punch file I/O
!  (2 ) error_mod.f    : Module containing NaN and other error check routines
!  (3 ) logical_mod.f  : Module containing GEOS-CHEM logical switches
!  (4 ) tracer_mod.f   : Module containing GEOS-CHEM tracer array STT etc.
!  (5 ) tracerid_mod.f : Module containing pointers to tracers & emissions
!  (6 ) pressure_mod.f : Module containing routines to compute P(I,J,L)
!
!  NOTES: 
!  (1 ) Routine "upbdflx_noy" now correctly reprocessed P(NOy) files from
!        /data/ctm/GEOS_4x5/pnoy_200106 or /data/ctm/GEOS_2x2.5/pnoy_200106.
!        (mje, bmy, 6/28/01)
!  (2 ) Updated comments (bmy, 9/4/01)
!  (3 ) Fixes for reading binary punch files of global size (bmy, 9/27/01)
!  (4 ) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (5 ) Removed obsolete commented out code from 7/01 (bmy, 11/26/01)
!  (6 ) Updated comments (bmy, 5/28/02)
!  (7 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in ordr
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (8 ) Now references "pressure_mod.f" (dsa, bdf, bmy, 8/21/02)
!  (9 ) Now references BXHEIGHT from "dao_mod.f".  Also deleted obsolete
!        code from 8/02.  Now references IDTNOx, IDTOX, from "tracerid_mod.f" 
!        instead of from "comtrid.h". (bmy, 11/6/02)
!  (10) Added driver routine DO_UPBDFLX.  Also added lat limits for 1x1 in
!        UPBDFLX_O3. (bmy, 3/14/03)
!  (11) Now references AD from "dao_mod.f" in UPBDFLX_NOY (bnd, bmy, 4/14/03)
!  (12) Added printout of O3 in Tg/yr in UPBDFLX_O3 (mje, bmy, 8/15/03)
!  (13) Change O3 flux for GEOS-3 to 500 Tg/yr in UPBDFLX_O3 (bmy, 9/15/03)
!  (14) Now references "tagged_ox_mod.f" (bmy, 8/19/03)
!  (15) Now activated parallel DO loops (bmy, 4/15/04)
!  (16) Now made IORD, JORD, KORD module variables.  Now added routine
!        SET_UPBDFLX.  Now added routine SET_TRANSPORT (bmy, 7/20/04)
!  (17) Bug fix for COMPAQ compiler.  Now supports 1x125 grid. (bmy, 12/1/04)
!******************************************************************************
!      
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "upbdflx_mod.f"
      !=================================================================
      
      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC  :: DO_UPBDFLX_ADJ
      PUBLIC  :: UPBDFLX_NOY_ADJ

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      INTEGER :: IORD, JORD, KORD

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
      
      SUBROUTINE DO_UPBDFLX_ADJ
!
!******************************************************************************
!  Subroutine DO_UPBDFLX is the driver routine for the stratospheric (upper-
!  boundary) routines for Ox and NOy. (bmy, 3/11/03, 7/20/04)
!  
!  NOTES:
!  (1 ) Removed IORD, JORD, KORD from the arg list.  Now references LPRT
!        from "logical_mod.f".  Now references ITS_A_FULLCHEM_SIM and
!        ITS_A_TAGOX_SIM from "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Add UPBDFLX_NOY_ADJ (dkh, 02/22/10) 
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE LOGICAL_MOD, ONLY : LPRT
      USE TRACER_MOD,  ONLY : ITS_A_FULLCHEM_SIM, ITS_A_TAGOX_SIM
      USE LINOZ_ADJ_MOD

#     include "CMN_SIZE"  ! Size parameters

      !=================================================================
      ! DO_UPBDFLX begins here!
      !=================================================================

      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         !---------------
         ! Fullchem run
         !---------------

         ! NOy from strat
         CALL UPBDFLX_NOY_ADJ( 1 )


         ! Ox from strat 
!         CALL UPBDFLX_O3
         !dbj changed for linoz
         CALL DO_LINOZ_ADJ

      ELSE IF ( ITS_A_TAGOX_SIM() ) THEN

         !---------------
         ! Tagged Ox run
         !---------------

         ! Ox from strat
!         CALL UPBDFLX_O3
         !dbj changed for linoz
         CALL DO_LINOZ_ADJ

      ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_UPBDFLX: after strat fluxes' )

      ! Return to calling program
      END SUBROUTINE DO_UPBDFLX_ADJ
!------------------------------------------------------------------------------

      SUBROUTINE UPBDFLX_NOY_ADJ( IFLAG )
!
!******************************************************************************
!  Subroutine UPBDFLX_NOY_ADJ is the adjoint of UPBDFLX_NOY. (dkh, 02/22/10) 
!
!  Based on forward model routine by  (qli, rvm, mje, bmy, 12/22/99, 8/4/06)
!
!  Arguments as input:
!  ===========================================================================
!  (1) IFLAG : IFLAG=1 will partition    [NOy] before transport
!              IFLAG=2 will re-partition [NOy] after  transport
!
!  NOTES:
!  (1 ) See forward model. 
!
!******************************************************************************
!      
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,     READ_BPCH2
      USE DAO_MOD,        ONLY : AD
      USE DIRECTORY_MOD,  ONLY : DATA_DIR
      USE ERROR_MOD,      ONLY : ERROR_STOP    
      USE TRACERID_MOD,   ONLY : IDTNOX,       IDTHNO3
      USE TIME_MOD,       ONLY : GET_TS_DYN,   GET_MONTH
      USE TIME_MOD,       ONLY : ITS_A_NEW_MONTH
      USE TIME_MOD,       ONLY : GET_NYMD
      USE TIME_MOD,       ONLY : GET_NYMDb
      USE TRACER_MOD,     ONLY : STT,          XNUMOLAIR
      USE TRANSFER_MOD,   ONLY : TRANSFER_ZONAL
      USE TROPOPAUSE_MOD, ONLY : GET_MIN_TPAUSE_LEVEL
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: IFLAG

      ! Local variables
      INTEGER              :: I, J, L, LMIN
      INTEGER, SAVE        :: LASTMONTH = -99

      TYPE (XPLEX)               :: DTDYN, AIRDENS, PNOY
      TYPE (XPLEX)               :: ARRAY(1,JGLOB,LGLOB)

      TYPE (XPLEX)               :: PNOY_ADJ
      INTEGER              :: MONTH_PRIOR

      ! Ratio of ( [NO] + [NO2] ) / [NOy] 
      TYPE (XPLEX), SAVE         :: XRATIO(JJPAR,LLPAR) 

      ! Arrays for P(NOY), NO, NO2, and HNO3 concentrations
      TYPE (XPLEX), SAVE         :: STRATPNOY(JJPAR,LLPAR)
      TYPE (XPLEX), SAVE         :: STRATNO(JJPAR,LLPAR)
      TYPE (XPLEX), SAVE         :: STRATNO2(JJPAR,LLPAR)
      TYPE (XPLEX), SAVE         :: STRATHNO3(JJPAR,LLPAR)

      ! For P(NOy) above 10 mb
      TYPE (XPLEX), SAVE         :: SPNOY10mb(JJPAR)

      ! TAU values for indexing the punch file 
      TYPE (XPLEX)               :: XTAU

      ! File Names
      CHARACTER (LEN=255)  :: FILENAME
      CHARACTER (LEN=255)  :: FILENAME2

      ! External functions
      TYPE (XPLEX),  EXTERNAL    :: BOXVL

      !=================================================================
      ! UPBDFLX_NOY_ADJ begins here! 
      !=================================================================

      ! Dynamic timestep [s]
      DTDYN = GET_TS_DYN() * 60d0

      !=================================================================
      ! IFLAG = 1: Before transport
      !
      ! If we have entered into a new month, read P(NOy), HNO3,
      ! NO, and NO2 from disk (binary punch file format).
      !=================================================================
      IF ( IFLAG == 1 ) THEN

         ! fwd:
         !IF ( ITS_A_NEW_MONTH() ) THEN
         ! adj:
         IF ( ITS_A_NEW_MONTH() .and. GET_NYMD() .ne. GET_NYMDb() ) THEN

            ! adj: calculate month prior 
            MONTH_PRIOR = GET_MONTH() - 1
            IF ( MONTH_PRIOR == 0 ) MONTH_PRIOR = 12

            ! fwd:
            ! TAU value corresponding to the beginning of this month
            !XTAU = GET_TAU0( GET_MONTH(), 1, 1985 ) 
            ! adj:
            ! TAU value corresponding to the beginning of previous month
            XTAU = GET_TAU0( MONTH_PRIOR, 1, 1985 ) 

            ! File containing P(NOy), NOx, HNO3 concentrations
            ! Now read corrected file from pnoy_200106/ subdir (bmy, 6/28/01)
            FILENAME  = TRIM( DATA_DIR )             // 
     &                  'pnoy_200106/pnoy_nox_hno3.' //
     &                  GET_NAME_EXT()   // '.'      // GET_RES_EXT()

            ! Echo filename to stdout
            WRITE( 6, 100 ) TRIM( FILENAME )
 100        FORMAT( '     - UPBDFLX_NOY: Reading ', a )
            
            ! P(NOy) in [v/v/s] is stored as tracer #1
            CALL READ_BPCH2( FILENAME, 'PNOY-L=$', 1,     
     &                       XTAU,      1,         JGLOB,     
     &                       LGLOB,     ARRAY,     QUIET=.TRUE. )

            ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize to (JJPAR,LLPAR)
            CALL TRANSFER_ZONAL( ARRAY(1,:,:), STRATPNOY )

            ! [HNO3] in [v/v] is stored as tracer #2
            CALL READ_BPCH2( FILENAME, 'PNOY-L=$', 2,     
     &                       XTAU,      1,         JGLOB,     
     &                       LGLOB,     ARRAY,     QUIET=.TRUE. )

            ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize to (JJPAR,LLPAR)
            CALL TRANSFER_ZONAL( ARRAY(1,:,:), STRATHNO3 )
         
            ! [NO] in [v/v] is stored as tracer #4
            CALL READ_BPCH2( FILENAME, 'PNOY-L=$', 4,     
     &                       XTAU,      1,         JGLOB,     
     &                       LGLOB,     ARRAY,     QUIET=.TRUE. )

            ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize to (JJPAR,LLPAR)
            CALL TRANSFER_ZONAL( ARRAY(1,:,:), STRATNO )

            ! [NO2] in [v/v] is stored as tracer #5
            CALL READ_BPCH2( FILENAME, 'PNOY-L=$', 5,     
     &                       XTAU,      1,         JGLOB,     
     &                       LGLOB,     ARRAY,     QUIET=.TRUE. )

            ! Cast from TYPE (XPLEX) to COMPLEX*16 and resize to (JJPAR,LLPAR)
            CALL TRANSFER_ZONAL( ARRAY(1,:,:), STRATNO2 )

            !===========================================================
            ! XRATIO is the ratio ( [NO] + [NO2] ) / [NOy],
            ! which is needed for the partitioning.
            ! XRATIO will be the same for a given month
            !===========================================================
            DO L = 1, LLPAR
            DO J = 1, JJPAR
               XRATIO(J,L) = ( STRATNO(J,L) + STRATNO2(J,L) ) /
     &                       ( STRATNO(J,L) + STRATNO2(J,L) + 
     &                         STRATHNO3(J,L) ) 
            ENDDO
            ENDDO
         ENDIF

         !==============================================================
         ! Initial partitioning of [NOy] to [NOx] and [HNO3], before 
         ! transport
         ! 
         ! We use zonal mean values for stratospheric P(NOy), [NO], 
         ! [NO2], and [HNO3] taken from Dylan Jones' & Hans Schneider's 
         ! 2-D model.  
         !
         ! Since P(NOy) above 10mb accounts for almost 50% of the total 
         ! stratospheric production, we also dump P(NOy) above 10 mb 
         ! into the top layer of the model.  These values are also 
         ! supplied to us by Dylan Jones.
         !
         ! We make the following assumptions:
         !
         !    (1) [NOx] = [NO] + [NO2] 
         !    (2) [NOy] = [NO] + [NO2] + [HNO3] = [NOx] + [HNO3]  
         !
         ! Therefore, in order to obtain [NOx] and [HNO3] from [NOy], 
         ! we must do the partitioning as follows:
         !
         !    (1) [NOy]   = P(NOy) + [NOx] + [HNO3]  
         !                = Production of NOy plus current 
         !                  concentrations of NOx and HNO3 in the 
         !                  given grid box
         !
         !    (2) XRATIO  = ( [NO] + [NO2] ) / [NOy]
         !
         !    (3) P(NOx)  = P(NOy) * XRATIO
         !  
         !    (4) P(HNO3) = P(NOy) * ( 1 - XRATIO ) 
         !
         ! XRATIO = ( [NO] + [NO2] ) / [NOy] approximates the true 
         ! ratio of [NOx] / [NOy], but is itself not the true ratio, 
         ! since our formulation of [NOy] neglects some additional 
         ! species (e.g. PAN, HNO4, N2O5, R4N2, PPN, PMN).
         !
         ! At some future point we may take the additional constituents 
         ! of [NOy] into account.  For now we proceed as outlined above.
         !==============================================================

         ! Minimum value of LPAUSE
         LMIN = GET_MIN_TPAUSE_LEVEL()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, PNOY )
!$OMP+PRIVATE( PNOY_ADJ      )
!$OMP+SCHEDULE( DYNAMIC )
         DO L = LMIN, LLPAR
         DO J = 1,    JJPAR
         DO I = 1,    IIPAR

            ! Skip over tropospheric boxes
            IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN 

               ! fwd code:
               !STT(I,J,L,IDTHNO3) = PNOY *
               !                     MAX( ( 1d0 -  XRATIO(J,L) ), 1d-20 )
               ! adj code:
               PNOY_ADJ = STT_ADJ(I,J,L,IDTHNO3) *
     &                              MAX( ( 1d0 -  XRATIO(J,L) ), 1d-20 )

             
               ! fwd code:
               !STT(I,J,L,IDTNOX)  = PNOY * XRATIO(J,L) 
               ! adj code: 
               PNOY_ADJ = PNOY_ADJ + XRATIO(J,L) * STT_ADJ(I,J,L,IDTNOX)

               ! fwd code:
               !PNOY = PNOY + STT(I,J,L,IDTNOX) + STT(I,J,L,IDTHNO3)
               ! adj code: note that STT gets overwritten in fwd code
               ! so STT_ADJ is not additive here. 
               STT_ADJ(I,J,L,IDTNOX)  = PNOY_ADJ
               STT_ADJ(I,J,L,IDTHNO3) = PNOY_ADJ



            ENDIF
         ENDDO   
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! IFLAG = 2: After transport
      !
      ! Repartition [NOy] after transport into [NOx] + [HNO3]
      !
      ! This repartitioning is necessary to avoid performing chemistry
      ! between the [NO2] and [HNO3] species.
      !
      ! The concentrations [NOx] and [HNO3] will have changed due to 
      ! transport, but the ratio used to partition them will be the 
      ! same.
      !=================================================================
      ELSE IF ( IFLAG == 2 ) THEN

         ! Minimum value of LPAUSE
         LMIN = GET_MIN_TPAUSE_LEVEL()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, PNOY )
!$OMP+PRIVATE( PNOY_ADJ      )
!$OMP+SCHEDULE( DYNAMIC )
         DO L = LMIN, LLPAR
         DO J = 1,    JJPAR
         DO I = 1,    IIPAR

            ! Skip over tropospheric boxes
            IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN

               ! Partition total [NOy] to [HNO3], units are [v/v]
               !STT(I,J,L,IDTHNO3) = PNOY *
               !                     MAX( ( 1d0 -  XRATIO(J,L) ), 1d-20 )
               PNOY_ADJ  = STT_ADJ(I,J,L,IDTHNO3) *
     &                              MAX( ( 1d0 -  XRATIO(J,L) ), 1d-20 )
 
               ! fwd code:
               !STT(I,J,L,IDTNOX)  = PNOY * XRATIO(J,L) 
               ! adj code:
               PNOY_ADJ = PNOY_ADJ + XRATIO(J,L) * STT_ADJ(I,J,L,IDTNOX)  

               ! fwd code:
               !PNOY = STT(I,J,L,IDTNOX) + STT(I,J,L,IDTHNO3)
               ! adj code:
               STT_ADJ(I,J,L,IDTNOX)  = PNOY_ADJ 
               STT_ADJ(I,J,L,IDTHNO3) = PNOY_ADJ
                

            ENDIF
         ENDDO   
         ENDDO
         ENDDO   
!$OMP END PARALLEL DO      

      ELSE 

         ! If IFLAG /= 1 or IFLAG /= 2, print an error message and stop
         CALL ERROR_STOP( 'IFLAG must be 1 or 2!',
     &                    'UPBDFLX_NOY_ADJ (upbdflx_adj_mod.f)' )

      ENDIF

      ! Return to calling program
      END SUBROUTINE UPBDFLX_NOY_ADJ

!------------------------------------------------------------------------------

      ! End of module
      END MODULE UPBDFLX_ADJ_MOD
