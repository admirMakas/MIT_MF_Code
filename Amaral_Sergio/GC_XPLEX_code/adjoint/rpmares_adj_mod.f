! $Id: rpmares_adj_mod.f,v 1.1 2009/09/09 06:12:55 daven Exp $
      MODULE RPMARES_ADJ_MOD
!
!******************************************************************************
!  Module RPMARES_ADJ_MOD is used to call the aerosol thermo adjoint 
!  subroutines (dkh, 09/08/09) 
! 
!  Module Routines:
!  ============================================================================
!  (1 ) DO_RPMARES_ADJ     : Driver which calls adjoint thermo routines
!  (2 ) adactcof           : adjoint of actcof 
!  (3 ) adawater           : adjoint of awater  
!  (4 ) adrpmares_11       : adjoint of rpmares exit=11
!  (5 ) adrpmares_12       : adjoint of rpmares exit=12
!  (6 ) adrpmares_2        : adjoint of rpmares exit=2
!  (7 ) adrpmares_3        : adjoint of rpmares exit=3
!  (8 ) adrpmares_4        : adjoint of rpmares exit=4
!  (9 ) adrpmares_6        : adjoint of rpmares exit=6 (old)
!  (10) adrpmares_7        : adjoint of rpmares exit=7
!  (11) adrpmares_8        : adjoint of rpmares exit=8
!  (12) adcubic            : adjoint of cubic 
!  (13) adrpmares_6_D5     : adjoint of rpmares exit=6 (correct)
!
!  GEOS-CHEM modules referenced by chemistry_mod.f
!  ============================================================================
!  (1 ) checkpt_mod        : Module w/ routines for checkpointing 
!  (2 ) dao_mod            : Module containing arrays for DAO met fields
!  (3 ) rpmares_mod        : Module w/ routines for aerosol thermodynamics 
!  (4 ) tracerid_mod       : Module containing pointers to tracers & emissions
!
!  NOTES:
!  (1 )  Updated to GCv8 (dkh, 09/09/09) 
!  
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
      
      SUBROUTINE DO_RPMARES_ADJ ( )
!
!******************************************************************************
!  Subroutine DO_RPMARES_ADJ is the driver for the TAMC generated adjoint of the
!  aerosol thermodynamic routine RPMARES 
!  (dkh, 8/27/04, 09/08/09) 
!
!  Passed via checkpoint_mod.f
!  ============================================================================
!
!  NOTES:
!  (1 ) Add error checking.  Add ADJ_NAN, FIRST, ADJ_NAN_COUNT, ADJ_MAX...
!       Check ADJ_STT for NAN and for large increases after calls to adrpmares.
!       Now reference IT_IS_NAN. (dkh, 02/08/05)
!  (2 ) Move paramters NCTOT and NPAR from CMN_ADJ to here. Replace many uses 
!        of NADJ with NCTOT.  Initialize ADJ_TMP after initializing ADJ_STT_LOCAL.
!        ADJ_STT_LOCAL is now dim = 8, not dim = NOBS, as NOBS may be much larger. 
!        Remove IS_DURING_OBSERVATION argument 
!        No longer force adjoints in this routine, do it in ???. (dkh, 03/03/05)        
!  (3 )  Now the ITS_TIME_FOR_CHEM section uses ADJ_STT [kg], so switch to
!         [ug/m3] for this portion, and switch back at the end. (dkh, 03/10/05)
!  (4 ) Replace the adjoint code for the high ratio case ( NRETURN = 6 ) with 
!        improved code that is more accurate and requires less checkpointing. 
!        (dkh, 06/01/05)
!  (5 ) Now reference ADJ_CONVERT_UNITS from dao_mod.f (dkh, 11/03/05)  
!  (6 ) Udpated to GCv8, renamed from ADJ_AEROSOL to DO_RPMARES_ADJ. 
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD, LFD
      USE ADJ_ARRAYS_MOD,  ONLY : STT_ADJ
      USE CHECKPT_MOD, 	   ONLY : RP_IN, RP_OUT
      USE DAO_MOD,         ONLY : AIRVOL
      USE ERROR_MOD,       ONLY : ERROR_STOP, IT_IS_NAN
      USE LOGICAL_ADJ_MOD, ONLY : LPRINTFD
      USE TRACERID_MOD,    ONLY : IDTSO4, IDTNH3, IDTNH4, IDTHNO3 
      USE TRACERID_MOD,    ONLY : IDTNIT
      USE TROPOPAUSE_MOD,  ONLY : ITS_IN_THE_STRAT


#     include "CMN_SIZE"    ! Size params

      ! Parameters
      INTEGER, PARAMETER  ::  MAX_ALLOWED_NAN      = 10
      INTEGER, PARAMETER  ::  MAX_ALLOWED_EXPLD    = 10
      TYPE (XPLEX),  PARAMETER::MAX_ALLOWED_INCREASE=xplex(10.0D10,0d0)

      INTEGER, PARAMETER  ::  NPAR  = 2  ! Number of input variables to RPMARES that are
                                         !  parameters, i.e. temp and rh

      INTEGER, PARAMETER  ::  NCTOT = 5  ! Number of input variables to RPMARES that are 
                                         !  total concentrations, = NRPIN - NPAR 

      ! Local variables
      TYPE (XPLEX)  ::  CTOT_P(NCTOT)       ! Same size as argument of the ad_rpmares routines
      TYPE (XPLEX)  ::  PAR_P(NPAR) 
      TYPE (XPLEX)  ::  ADJ_STT_LOCAL(8)    ! Same size as argument of the ad_rpmares routines
      TYPE (XPLEX)  ::  ADJ_CTOT(NCTOT)  
      INTEGER ::  I, J, L, N
      INTEGER ::  NRETURN
      TYPE (XPLEX)  ::  ADJ_TMP(NCTOT)     ! Temp storage for resetting bad adjs to original value
      TYPE (XPLEX)  ::  MAX_ADJ_TMP        ! Temp max value used for error checking
      LOGICAL ::  ADJ_NAN = .FALSE. 
      INTEGER ::  ADJ_NAN_COUNT, ADJ_EXPLD_COUNT 
      LOGICAL, SAVE :: FIRST = .TRUE.
      TYPE (XPLEX)  ::  AVOL 
     
      !================================================================
      ! DO_RPMARES_ADJ begins here!
      !================================================================
   
      ! Initialize ADJ_NAN_COUNT the first time through
      IF ( FIRST ) THEN
         ADJ_NAN_COUNT   = 0 
         ADJ_EXPLD_COUNT = 0 
      
         FIRST         = .FALSE.
      ENDIF

      ! Save maximum adjoint for error checking later
      MAX_ADJ_TMP = MAXVAL( ABS(STT_ADJ) )
      
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,     L,       N   )
!$OMP+PRIVATE( CTOT_P,  PAR_P, NRETURN, ADJ_STT_LOCAL )
!$OMP+PRIVATE( ADJ_TMP, ADJ_CTOT  )
!$OMP+PRIVATE( AVOL )
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip if we are in the stratosphere (bmy, 4/3/08)
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE

         ! air volume 
         AVOL  = AIRVOL(I,J,L)

         ! Load IN from RP_IN 
         CTOT_P(:) = RP_IN(I,J,L,1:NCTOT) 

         ! Load parameters from RP_IN 
         PAR_P(:)  = RP_IN(I,J,L,6:7) 
 
         ! Find out where RPMARES exited during the forward run
         NRETURN   = RP_OUT(I,J,L,9)


         ! Initialize the independent  and dependent variables to 0
         ADJ_CTOT(:)      = 0.D0
         ADJ_STT_LOCAL(:) = 0.D0

         ! Copy current value of ADJ variable to ADJ_STT_LOCAL 
         ! Always update the local adjoint input to the current adjoint tracer 
         ! values
         ADJ_STT_LOCAL(3) = STT_ADJ(I,J,L,IDTNIT) * AVOL * 1.d-9 
         ADJ_STT_LOCAL(5) = STT_ADJ(I,J,L,IDTNH4) * AVOL * 1.d-9  
         ADJ_STT_LOCAL(7) = STT_ADJ(I,J,L,IDTHNO3) * AVOL * 1.d-9
         ADJ_STT_LOCAL(8) = STT_ADJ(I,J,L,IDTNH3) * AVOL * 1.d-9
         ! Since thermo doesn't modify total sulfate, don't need to 
         ! pass it initial adjoint values for SO4 ?
         ADJ_STT_LOCAL(6) = 0d0 

         ! The forcing for these species is also zero
         ADJ_STT_LOCAL(1) = 0.d0
         ADJ_STT_LOCAL(2) = 0.d0
         ADJ_STT_LOCAL(4) = 0.d0

         ! Store original values in ADJ_TMP
         ADJ_TMP(1) = STT_ADJ(I,J,L,IDTSO4)  
         ADJ_TMP(2) = STT_ADJ(I,J,L,IDTHNO3)
         ADJ_TMP(3) = STT_ADJ(I,J,L,IDTNH3)
         ADJ_TMP(4) = STT_ADJ(I,J,L,IDTNIT)
         ADJ_TMP(5) = STT_ADJ(I,J,L,IDTNH4)

         IF ( LPRINTFD 
     &      .and. J == JFD .AND. L == LFD .AND. I == IFD) THEN 
            print*, 'before ', CTOT_P, PAR_P, ADJ_CTOT, ADJ_STT_LOCAL
            print*, 'NRETURN = ', NRETURN
         ENDIF 
         
         !============================================================
         ! CALCULATE ADJOINT THERMO 
         !
         ! The thermodynamic routine is broken into several regimes. 
         ! The regime from the forward calculation is marked by the 
         ! NRETURN flag.  Use this flag to push the adjoint calculation
         ! into the same regime.  
         !============================================================         

         IF (NRETURN == 1) THEN
            ! the adjoint variables are unchanged in this case
            ADJ_CTOT(:) = ADJ_TMP(:)

         ELSEIF (NRETURN == 2) THEN
            CALL adrpmares_2( CTOT_P, PAR_P, ADJ_CTOT, ADJ_STT_LOCAL,
     &                     I,      J,   L )
         ELSEIF (NRETURN == 3 .OR. NRETURN == 5) THEN
            CALL adrpmares_3( CTOT_P, PAR_P, ADJ_CTOT, ADJ_STT_LOCAL,
     &                     I,      J,   L )
         ELSEIF (NRETURN == 4) THEN
            CALL adrpmares_4( CTOT_P, PAR_P, ADJ_CTOT, ADJ_STT_LOCAL,
     &                     I,      J,   L )
         ELSEIF (NRETURN == 6) THEN
            CALL adrpmares_6_D5( CTOT_P, PAR_P, ADJ_CTOT, ADJ_STT_LOCAL,
     &                     I,      J,   L )
         ELSEIF (NRETURN == 7) THEN
            CALL adrpmares_7( CTOT_P, PAR_P, ADJ_CTOT, ADJ_STT_LOCAL,
     &                     I,      J,   L )
         ELSEIF (NRETURN == 8 .OR. NRETURN == 9
     &      .OR. NRETURN == 10) THEN
            CALL adrpmares_8( CTOT_P, PAR_P, ADJ_CTOT, ADJ_STT_LOCAL,
     &                      I,      J,   L )
         ELSEIF (NRETURN == 11) THEN
            CALL adrpmares_11( CTOT_P, PAR_P, ADJ_CTOT, ADJ_STT_LOCAL,
     &                      I,      J,   L )
         ELSEIF (NRETURN == 12) THEN
            CALL adrpmares_12( CTOT_P, PAR_P, ADJ_CTOT, ADJ_STT_LOCAL,
     &                      I,      J,   L )
         ELSE
            print*, ' NRETURN = ', NRETURN , I, J, L 
            CALL ERROR_STOP
     &                 ('ERROR: NRETURN ill defined ','ADJ_AEROSOL')
         ENDIF

         IF ( LPRINTFD .AND. 
     &        J == JFD .AND. L == LFD .AND. I == IFD) THEN 
            print*, 'after ', CTOT_P, PAR_P, ADJ_CTOT, ADJ_STT_LOCAL
         ENDIF 

         ! Check for NAN. 
         DO N = 1, NCTOT

            IF ( IT_IS_NAN( ADJ_CTOT(N) ) ) THEN

               ! Echo location of NAN (probably leave this commented out
               ! unless you are getting lots of ADJ_NAN warnings
               !WRITE(6,*) 'FOUND A NAN AT I,J,L,N = ',I,J,L,N

!$OMP CRITICAL
               ! Set ADJ_NAN flag so that a warning is echod to screen
               ADJ_NAN = .TRUE.
!$OMP END CRITICAL

               ! Replace the NAN with the original value and continue
               ADJ_CTOT(N) = ADJ_TMP(N)
            ENDIF

         ENDDO

         ! Update ADJ_STT array
         STT_ADJ(I,J,L,IDTHNO3) = ADJ_CTOT(2) * 1.d9 / AVOL 
         STT_ADJ(I,J,L,IDTNH3)  = ADJ_CTOT(3) * 1.d9 / AVOL
         STT_ADJ(I,J,L,IDTNIT)  = ADJ_CTOT(4) * 1.d9 / AVOL
         STT_ADJ(I,J,L,IDTNH4)  = ADJ_CTOT(5) * 1.d9 / AVOL
         
         ! Becuase we don't initiate the sulfate adjoint with STT_ADJ,
         ! do not overwrite. 
         STT_ADJ(I,J,L,IDTSO4)  = STT_ADJ(I,J,L,IDTSO4) 
     &                          + ADJ_CTOT(1) * 1.d9 / AVOL

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO


      ! Error checking
      IF ( ADJ_NAN ) THEN
     
         ! Echo a warning to the screen
         WRITE(6,*)
     &             ' *** - WARNING: ADJ_NAN in routine ADJ_AEROSOL'

         ! keep track of how many times NANs have occured
         ADJ_NAN_COUNT = ADJ_NAN_COUNT + 1 

         IF ( ADJ_NAN_COUNT > MAX_ALLOWED_NAN )
     &       CALL ERROR_STOP('Too many NANs', 'ADJ_AEROSOL')

      ENDIF

      ! More error checking: warn of exploding adjoit values, except
      ! the first jump up from zero (MAX_ADJ_TMP = 0 first few times)
      IF ( MAXVAL(ABS(STT_ADJ)) > (MAX_ADJ_TMP * MAX_ALLOWED_INCREASE) 
     &   .AND. ( MAX_ADJ_TMP > 0d0 )  ) THEN

         WRITE(6,*)' *** - WARNING: EXPLODING adjoints in ADJ_AEROSOL' 
         WRITE(6,*)' *** - MAX(ADJ_STT) before = ',MAX_ADJ_TMP 
         WRITE(6,*)' *** - MAX(ADJ_STT) after  = ',MAXVAL(ABS(STT_ADJ)) 
     
         ADJ_EXPLD_COUNT = ADJ_EXPLD_COUNT + 1
         
         IF (ADJ_EXPLD_COUNT > MAX_ALLOWED_EXPLD )
     &      CALL ERROR_STOP('Too many exploding adjoints',
     &                       'ADJ_AEROSOL, adjoint_mod.f')

       ENDIF 

      ! Return to calling program 
      END SUBROUTINE DO_RPMARES_ADJ
!------------------------------------------------------------------------------
   
C                           DISCLAIMER
C
C   This file was generated by TAMC version 5.3.2
C
C   THE AUTHOR DOES NOT MAKE  ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR  RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
C   OR USEFULNESS  OF ANY INFORMATION OR PROCESS  DISCLOSED, OR REPRESENTS
C   THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
C   THIS CODE IS FOR NON-PROFIT-ORIENTATED ACADEMIC RESEARCH AND EDUCATION
C   ONLY.  ANY COMMERCIAL OR OTHER PROFIT-ORIENTATED USE OR  EVALUATION IS
C   STRICTLY  FORBIDDEN.  PASSING  THIS CODE  TO  ANY THIRD  PARTY IS  NOT
C   ALLOWED.
C
C   FOR COMMERCIAL OR  OTHER PROFIT-ORIENTATED APPLICATIONS PLEASE CONTACT
C   info@FastOpt.com
C
      subroutine adactcof( cat, an, adcat, adan, adgama )
C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================

C==============================================
C define parameters                            
C==============================================
      integer nan
      parameter ( nan = 3 )
      integer ncat
      parameter ( ncat = 2 )

C==============================================
C define common blocks
C==============================================
      common /actcofsv/ pname, xmsg, zm, zp
      character*(16)  :: pname = ' driver program name'
      character*(120)  :: xmsg = ' '
      TYPE (XPLEX)  :: zm(nan) = (/xplex(2.d0,0d0)
     & , xplex(1.d0,0d0),xplex(1.d0,0d0)/)
      TYPE (XPLEX)  :: zp(ncat) = (/xplex(1.d0,0d0),
     & xplex(1.d0,0d0)/)

C==============================================
C define arguments
C==============================================
      TYPE (XPLEX) adan(nan)
      TYPE (XPLEX) adcat(ncat)
      TYPE (XPLEX) adgama(ncat,nan)
      TYPE (XPLEX) an(nan)
      TYPE (XPLEX) cat(ncat)

C==============================================
C define local variables
C==============================================
      TYPE (XPLEX) adbgama(ncat,nan)
      TYPE (XPLEX) adf1(nan)
      TYPE (XPLEX) adf2(ncat)
      TYPE (XPLEX) adfgama
      TYPE (XPLEX) adi
      TYPE (XPLEX) adlgama0(ncat,nan)
      TYPE (XPLEX) adm(ncat,nan)
      TYPE (XPLEX) adsri
      TYPE (XPLEX) adta
      TYPE (XPLEX) adtc
      TYPE (XPLEX) adtexpv
      TYPE (XPLEX) adtrm
      TYPE (XPLEX) adtwoi
      TYPE (XPLEX) adtwosri
      TYPE (XPLEX) adx(ncat,nan)
      TYPE (XPLEX) ady(nan,ncat)
      TYPE (XPLEX) adzot1
      TYPE (XPLEX) beta0(ncat,nan)
      TYPE (XPLEX) beta1(ncat,nan)
      TYPE (XPLEX) bgama(ncat,nan)
      TYPE (XPLEX) cgama(ncat,nan)
      integer exit
      TYPE (XPLEX) f1(nan)
      TYPE (XPLEX) f2(ncat)
      TYPE (XPLEX) fgama
      TYPE (XPLEX) i
      integer ian
      integer icat
      integer ip1
      integer ip2
      TYPE (XPLEX) lgama0(ncat,nan)
      TYPE (XPLEX) m(ncat,nan)
      TYPE (XPLEX) sri
      TYPE (XPLEX) ta
      TYPE (XPLEX) tb
      TYPE (XPLEX) tc
      TYPE (XPLEX) texpv
      TYPE (XPLEX) trm
      TYPE (XPLEX) twoi
      TYPE (XPLEX) twosri
      TYPE (XPLEX) v1(ncat,nan)
      TYPE (XPLEX) v2(ncat,nan)
      TYPE (XPLEX) x(ncat,nan)
      TYPE (XPLEX) y(nan,ncat)
      TYPE (XPLEX) zbar
      TYPE (XPLEX) zbar2
      TYPE (XPLEX) zot1

C==============================================
C define data
C==============================================
      data beta0(1,1)/xplex(2.98d-2,0d0)/
      data beta1(1,1)/xplex(0.0d0,0d0)/
      data cgama(1,1)/xplex(4.38d-2,0d0)/
      data beta0(1,2)/xplex(1.2556d-1,0d0)/
      data beta1(1,2)/xplex(2.8778d-1,0d0)/
      data cgama(1,2)/xplex(-5.59d-3,0d0)/
      data beta0(1,3)/xplex(2.0651d-1,0d0)/
      data beta1(1,3)/xplex(5.556d-1,0d0)/
      data cgama(1,3)/xplex(0.0d0,0d0)/
      data beta0(2,1)/xplex(4.6465d-2,0d0)/
      data beta1(2,1)/xplex(-0.54196d0,0d0)/
      data cgama(2,1)/xplex(-1.2683d-3,0d0)/
      data beta0(2,2)/xplex(-7.26224d-3,0d0)/
      data beta1(2,2)/xplex(-1.168858d0,0d0)/
      data cgama(2,2)/xplex(3.51217d-5,0d0)/
      data beta0(2,3)/xplex(4.494d-2,0d0)/
      data beta1(2,3)/xplex(2.3594d-1,0d0)/
      data cgama(2,3)/xplex(-2.962d-3,0d0)/
      data v1(1,1),v2(1,1)/xplex(2.0d0,0d0),xplex(1.0d0,0d0)/
      data v1(2,1),v2(2,1)/xplex(2.0d0,0d0),xplex(1.0d0,0d0)/
      data v1(1,2),v2(1,2)/xplex(1.0d0,0d0),xplex(1.0d0,0d0)/
      data v1(2,2),v2(2,2)/xplex(1.0d0,0d0),xplex(1.0d0,0d0)/
      data v1(1,3),v2(1,3)/xplex(1.0d0,0d0),xplex(1.0d0,0d0)/
      data v1(2,3),v2(2,3)/xplex(1.0d0,0d0),xplex(1.0d0,0d0)/

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      do ip2 = 1, nan
        do ip1 = 1, ncat
          adbgama(ip1,ip2) = 0.
        end do
      end do
      do ip1 = 1, nan
        adf1(ip1) = 0.
      end do
      do ip1 = 1, ncat
        adf2(ip1) = 0.
      end do
      adfgama = 0.
      adi = 0.
      do ip2 = 1, nan
        do ip1 = 1, ncat
          adlgama0(ip1,ip2) = 0.
        end do
      end do
      do ip2 = 1, nan
        do ip1 = 1, ncat
          adm(ip1,ip2) = 0.
        end do
      end do
      adsri = 0.
      adta = 0.
      adtc = 0.
      adtexpv = 0.
      adtrm = 0.
      adtwoi = 0.
      adtwosri = 0.
      do ip2 = 1, nan
        do ip1 = 1, ncat
          adx(ip1,ip2) = 0.
        end do
      end do
      do ip2 = 1, ncat
        do ip1 = 1, nan
          ady(ip1,ip2) = 0.
        end do
      end do
      adzot1 = 0.

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      exit = 0
      if (exit .eq. 0) then
        i = 0.d0
      endif
      do icat = 1, ncat
        if (exit .eq. 0) then
          i = i+cat(icat)*zp(icat)*zp(icat)
        endif
      end do
      do ian = 1, nan
        if (exit .eq. 0) then
          i = i+an(ian)*zm(ian)*zm(ian)
        endif
      end do
      if (exit .eq. 0) then
        i = 0.5d0*i
      endif
      if (i .eq. 0.d0) then
        if (exit .eq. 0) then
          exit = 1
        endif
      endif
      if (exit .eq. 0) then
        sri = sqrt(i)
        twosri = 2.d0*sri
        twoi = 2.d0*i
        texpv = 1.d0-exp(-twosri)*(1.d0+twosri-twoi)
        zot1 = 0.511d0*sri/(1.d0+sri)
        fgama = -(0.392d0*(sri/(1.d0+1.2d0*sri)+2.d0/1.2d0*log(1.d0+
     $1.2d0*sri)))
        do icat = 1, ncat
          do ian = 1, nan
            bgama(icat,ian) = 2.d0*beta0(icat,ian)+2.d0*beta1(icat,ian)/
     $(4.d0*i)*texpv
            m(icat,ian) = (cat(icat)**v1(icat,ian)*an(ian)**v2(icat,ian)
     $)**(1.d0/(v1(icat,ian)+v2(icat,ian)))
            lgama0(icat,ian) = (zp(icat)*zm(ian)*fgama+m(icat,ian)*2.d0*
     $v1(icat,ian)*v2(icat,ian)/(v1(icat,ian)+v2(icat,ian))*bgama(icat,
     $ian)+m(icat,ian)*m(icat,ian)*2.d0*(v1(icat,ian)*v2(icat,ian))**
     $1.5d0/(v1(icat,ian)+v2(icat,ian))*cgama(icat,ian))/2.302585093d0
          end do
        end do
        do ian = 1, nan
          do icat = 1, ncat
            zbar = (zp(icat)+zm(ian))*0.5d0
            zbar2 = zbar*zbar
            y(ian,icat) = zbar2*an(ian)/i
            x(icat,ian) = zbar2*cat(icat)/i
          end do
        end do
        do ian = 1, nan
          f1(ian) = 0.d0
          do icat = 1, ncat
            f1(ian) = f1(ian)+x(icat,ian)*lgama0(icat,ian)+zot1*zp(icat)
     $*zm(ian)*x(icat,ian)
          end do
        end do
        do icat = 1, ncat
          f2(icat) = 0.d0
          do ian = 1, nan
            f2(icat) = f2(icat)+y(ian,icat)*lgama0(icat,ian)+zot1*
     $zp(icat)*zm(ian)*y(ian,icat)
          end do
        end do
        do ian = 1, nan
          adta = 0.
          adtc = 0.
          adtrm = 0.
          do icat = 1, ncat
            adta = 0.
            adtc = 0.
            adtrm = 0.
            ta = -(zot1*zp(icat)*zm(ian))
            tb = zp(icat)*zm(ian)/(zp(icat)+zm(ian))
            tc = f2(icat)/zp(icat)+f1(ian)/zm(ian)
            trm = ta+tb*tc
            if (trm .gt. 30.d0) then
              adgama(icat,ian) = 0.
            else
              adtrm = adtrm+adgama(icat,ian)*10.d0**trm*dlog(10.d0)
              adgama(icat,ian) = 0.
            endif
            adta = adta+adtrm
            adtc = adtc+adtrm*tb
            adtrm = 0.
            adf1(ian) = adf1(ian)+adtc/zm(ian)
            adf2(icat) = adf2(icat)+adtc/zp(icat)
            adtc = 0.
            adzot1 = adzot1-adta*zp(icat)*zm(ian)
            adta = 0.
          end do
        end do
        do icat = 1, ncat
          do ian = 1, nan
            adlgama0(icat,ian) = adlgama0(icat,ian)+adf2(icat)*y(ian,
     $icat)
            ady(ian,icat) = ady(ian,icat)+adf2(icat)*(lgama0(icat,ian)+
     $zot1*zp(icat)*zm(ian))
            adzot1 = adzot1+adf2(icat)*zp(icat)*zm(ian)*y(ian,icat)
          end do
          adf2(icat) = 0.
        end do
        do ian = 1, nan
          do icat = 1, ncat
            adlgama0(icat,ian) = adlgama0(icat,ian)+adf1(ian)*x(icat,
     $ian)
            adx(icat,ian) = adx(icat,ian)+adf1(ian)*(lgama0(icat,ian)+
     $zot1*zp(icat)*zm(ian))
            adzot1 = adzot1+adf1(ian)*zp(icat)*zm(ian)*x(icat,ian)
          end do
          adf1(ian) = 0.
        end do
        do ian = 1, nan
          do icat = 1, ncat
            zbar = (zp(icat)+zm(ian))*0.5d0
            zbar2 = zbar*zbar
            adcat(icat) = adcat(icat)+adx(icat,ian)*(zbar2/i)
            adi = adi-adx(icat,ian)*(zbar2*cat(icat)/(i*i))
            adx(icat,ian) = 0.
            adan(ian) = adan(ian)+ady(ian,icat)*(zbar2/i)
            adi = adi-ady(ian,icat)*(zbar2*an(ian)/(i*i))
            ady(ian,icat) = 0.
          end do
        end do
        do icat = 1, ncat
          do ian = 1, nan
            bgama(icat,ian) = 2.d0*beta0(icat,ian)+2.d0*beta1(icat,ian)/
     $(4.d0*i)*texpv
            m(icat,ian) = (cat(icat)**v1(icat,ian)*an(ian)**v2(icat,ian)
     $)**(1.d0/(v1(icat,ian)+v2(icat,ian)))
            adbgama(icat,ian) = adbgama(icat,ian)+adlgama0(icat,ian)*
     $(m(icat,ian)*(2.d0*v1(icat,ian)*v2(icat,ian)/(v1(icat,ian)+
     $v2(icat,ian)))/2.302585093d0)
            adfgama = adfgama+adlgama0(icat,ian)*(zp(icat)*zm(ian)/
     $2.302585093d0)
            adm(icat,ian) = adm(icat,ian)+adlgama0(icat,ian)*((2.d0*
     $v1(icat,ian)*v2(icat,ian)/(v1(icat,ian)+v2(icat,ian))*bgama(icat,
     $ian)+2*m(icat,ian)*2.d0*(v1(icat,ian)*v2(icat,ian))**1.5d0/
     $(v1(icat,ian)+v2(icat,ian))*cgama(icat,ian))/2.302585093d0)
            adlgama0(icat,ian) = 0.

            ! The next two IF statements added to avoid divide by zero
            ! segmentation fault (dkh)
            IF (adm(icat,ian)*cat(icat)**v1(icat,ian)*
     $v2(icat,ian)*an(ian) .NE. 0.D0 .AND.
     $ (cat(icat)**v1(icat,ian)*an(ian)**v2(icat,ian)) .NE. 0.D0)
     $      adan(ian) = adan(ian)+adm(icat,ian)*cat(icat)**v1(icat,ian)*
     $v2(icat,ian)*an(ian)**(v2(icat,ian)-1)*1.d0/(v1(icat,ian)+v2(icat,
     $ian))*(cat(icat)**v1(icat,ian)*an(ian)**v2(icat,ian))**(1.d0/
     $(v1(icat,ian)+v2(icat,ian))-1)
            IF (adm(icat,ian)*v1(icat,ian)*
     $cat(icat) .NE. 0.D0 .AND. (cat(icat)**v1(icat,ian)*an(ian)**
     $v2(icat,ian)) .NE. 0.D0)
     $      adcat(icat) = adcat(icat)+adm(icat,ian)*v1(icat,ian)*
     $cat(icat)**(v1(icat,ian)-1)*an(ian)**v2(icat,ian)*1.d0/(v1(icat,
     $ian)+v2(icat,ian))*(cat(icat)**v1(icat,ian)*an(ian)**v2(icat,ian))
     $**(1.d0/(v1(icat,ian)+v2(icat,ian))-1)
            adm(icat,ian) = 0.
            adi = adi-adbgama(icat,ian)*8*beta1(icat,ian)/(16*i*i)*texpv
            adtexpv = adtexpv+adbgama(icat,ian)*(2.d0*beta1(icat,ian)/
     $(4.d0*i))
            adbgama(icat,ian) = 0.
          end do
        end do
        adsri = adsri-0.392d0*adfgama*(1/(1.d0+1.2d0*sri)-1.2d0*sri/
     $((1.d0+1.2d0*sri)*(1.d0+1.2d0*sri))+2*(1./(1.d0+1.2d0*sri)))
        adfgama = 0.
        adsri = adsri+adzot1*(0.511d0/(1.d0+sri)-0.511d0*sri/((1.d0+sri)
     $*(1.d0+sri)))
        adzot1 = 0.
        adtwoi = adtwoi+adtexpv*exp(-twosri)
        adtwosri = adtwosri-adtexpv*(exp(-twosri)-(1.d0+twosri-twoi)*
     $exp(-twosri))
        adtexpv = 0.
        adi = adi+2*adtwoi
        adtwoi = 0.
        adsri = adsri+2*adtwosri
        adtwosri = 0.
        adi = adi+adsri*(1./(2.*sqrt(i)))
        adsri = 0.
      endif
      exit = 0
      if (i .eq. 0.d0) then
        if (exit .eq. 0) then
          do ian = 1, nan
            do icat = 1, ncat
              adgama(icat,ian) = 0.
            end do
          end do
        endif
      endif
      if (exit .eq. 0) then
        adi = 0.5d0*adi
      endif
      do ian = 1, nan
        if (exit .eq. 0) then
          adan(ian) = adan(ian)+adi*zm(ian)*zm(ian)
        endif
      end do
      do icat = 1, ncat
        if (exit .eq. 0) then
          adcat(icat) = adcat(icat)+adi*zp(icat)*zp(icat)
        endif
      end do

      end SUBROUTINE ADACTCOF 
!-----------------------------------------------------------------------------

      subroutine adawater( irhx, mso4, mnh4, mno3, admso4, admnh4, 
     $admno3, adwh2o )
  
      ! Reference other f90 modules
      USE RPMARES_MOD, 	ONLY : POLY4, POLY6 

C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================

C==============================================
C define parameters                            
C==============================================
      TYPE (XPLEX) mwnh4
      parameter ( mwnh4 = xplex(18.0985d0,0d0) )
      TYPE (XPLEX) mwso4
      parameter ( mwso4 = xplex(96.0636d0,0d0) )
      TYPE (XPLEX) mw2
      parameter ( mw2 = xplex(mwso4%r+2.d0*mwnh4%r,0d0) )
      TYPE (XPLEX) mwno3
      parameter ( mwno3 = xplex(62.0649d0,0d0) )
      TYPE (XPLEX) mwano3
      parameter ( mwano3 =xplex( mwno3%r+mwnh4%r,0d0) )

C==============================================
C define common blocks
C==============================================
      common /awatersv/ c0, c1, c15, c2, kno3, kso4
      TYPE (XPLEX)  :: c0(4) = (/xplex(0.798079d0,0d0),
     & xplex(-1.574367d0,0d0),xplex(2.536686d0,0d0),
     & xplex(-1.735297d0,0d0)/)
      TYPE (XPLEX)  :: c1(4) = (/xplex(0.9995178d0,0d0),
     & xplex(-0.7952896d0,0d0),xplex(0.99683673d0,0d0),
     & xplex(-1.143874d0,0d0)/)
      TYPE (XPLEX)  :: c15(4) = (/xplex(1.697092d0,0d0),
     & xplex(-4.045936d0,0d0),xplex(5.833688d0,0d0),
     & xplex(-3.463783d0,0d0)/)
      TYPE (XPLEX)  :: c2(4) = (/xplex(2.085067d0,0d0),
     & xplex(-6.024139d0,0d0),xplex(8.967967d0,0d0),
     & xplex(-5.002934d0,0d0)/)
      TYPE (XPLEX)  ::kno3(6)=(/xplex(0.2906d0,0d0),
     & xplex(6.83665d0,0d0),xplex(-26.9093d0,0d0),
     $ xplex(46.6983d0,0d0),xplex(-38.803d0,0d0),
     & xplex(11.8837d0,0d0)/)
      TYPE (XPLEX)  ::kso4(6)=(/xplex(2.27515d0,0d0),
     & xplex(-11.147d0,0d0),xplex(36.3369d0,0d0),
     & xplex(-64.2134d0,0d0),xplex(56.8341d0,0d0),
     & xplex(-20.0953d0,0d0)/)

C==============================================
C define arguments
C==============================================
      TYPE (XPLEX) admnh4
      TYPE (XPLEX) admno3
      TYPE (XPLEX) admso4
      TYPE (XPLEX) adwh2o
      integer irhx
      TYPE (XPLEX) mnh4
      TYPE (XPLEX) mno3
      TYPE (XPLEX) mso4

C==============================================
C define local variables
C==============================================
      TYPE (XPLEX) adawc
      TYPE (XPLEX) adtnh4
      TYPE (XPLEX) adtno3
      TYPE (XPLEX) adtso4
      TYPE (XPLEX) adx
      TYPE (XPLEX) ady
      TYPE (XPLEX) ady40
      TYPE (XPLEX) adyc
      TYPE (XPLEX) aw
      TYPE (XPLEX) awc
      integer irh
      TYPE (XPLEX) mfs0
      TYPE (XPLEX) mfs1
      TYPE (XPLEX) mfs15
      TYPE (XPLEX) mfsno3
      TYPE (XPLEX) mfsso4
      TYPE (XPLEX) tnh4
      TYPE (XPLEX) tno3
      TYPE (XPLEX) tso4
      TYPE (XPLEX) x
      TYPE (XPLEX) y
      TYPE (XPLEX) y0
      TYPE (XPLEX) y1
      TYPE (XPLEX) y140
      TYPE (XPLEX) y15
      TYPE (XPLEX) y1540
      TYPE (XPLEX) y2
      TYPE (XPLEX) y3
      TYPE (XPLEX) y40
      TYPE (XPLEX) yc

!C==============================================
!C define external procedures and functions
!C==============================================
!      TYPE (XPLEX) poly4
!      external poly4
!      TYPE (XPLEX) poly6
!      external poly6
!
C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adawc = 0.
      adtnh4 = 0.
      adtno3 = 0.
      adtso4 = 0.
      adx = 0.
      ady = 0.
      ady40 = 0.
      adyc = 0.

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      irh = irhx
      irh = max(1,irh)
      irh = min(irh,100)
      aw = xplx(irh)/100.d0
      tso4 = max(mso4,0.d0)
      tnh4 = max(mnh4,0.d0)
      tno3 = max(mno3,0.d0)
      x = 0.d0
      if (tso4 .gt. 0.d0) then
        x = tnh4/tso4
      else
        if (tno3 .gt. 0.d0 .and. tnh4 .gt. 0.d0) then
          x = 10.d0
        endif
      endif
      if (x .lt. 1.d0) then
        mfs0 = poly4(c0,aw)
        mfs1 = poly4(c1,aw)
        y0 = (1.d0-mfs0)/mfs0
        y1 = (1.d0-mfs1)/mfs1
        y = (1.d0-x)*y0+x*y1
      else if (x .lt. 1.5d0) then
        if (irh .ge. 40) then
          mfs1 = poly4(c1,aw)
          mfs15 = poly4(c15,aw)
          y1 = (1.d0-mfs1)/mfs1
          y15 = (1.d0-mfs15)/mfs15
          y = 2.d0*(y1*(1.5d0-x)+y15*(x-1.d0))
        else
          awc = 0.8d0*(x-1.d0)
          y = 0.d0
          if (aw .ge. awc) then
            mfs1 = poly4(c1,xplx(0.4d0))
            mfs15 = poly4(c15,xplx(0.4d0))
            y140 = (1.d0-mfs1)/mfs1
            y1540 = (1.d0-mfs15)/mfs15
            y40 = 2.d0*(y140*(1.5d0-x)+y1540*(x-1.d0))
            yc = 2.d0*y1540*(x-1.d0)
            y = y40-(y40-yc)*(0.4d0-aw)/(0.4d0-awc)
          endif
        endif
      else if (x .lt. 2.d0) then
        y = 0.d0
        if (irh .ge. 40) then
          mfs15 = poly4(c15,aw)
          y15 = (1.d0-mfs15)/mfs15
          mfsso4 = poly6(kso4,aw)
          y2 = (1.d0-mfsso4)/mfsso4
          y = 2.d0*(y15*(2.d0-x)+y2*(x-1.5d0))
        endif
      else
        y2 = 0.d0
        y3 = 0.d0
        if (irh .ge. 40) then
          mfsso4 = poly6(kso4,aw)
          mfsno3 = poly6(kno3,aw)
          y2 = (1.d0-mfsso4)/mfsso4
          y3 = (1.d0-mfsno3)/mfsno3
        endif
      endif
      if (x .lt. 2.d0) then
        adtnh4 = adtnh4+adwh2o*y*mwnh4
        adtso4 = adtso4+adwh2o*y*mwso4
        ady = ady+adwh2o*(tso4*mwso4+mwnh4*tnh4)
        adwh2o = 0.
      else
        adtno3 = adtno3+adwh2o*y3*mwano3
        adtso4 = adtso4+adwh2o*y2*mw2
        adwh2o = 0.
      endif
      if (x .lt. 1.d0) then
        adx = adx+ady*((-y0)+y1)
        ady = 0.
      else if (x .lt. 1.5d0) then
        if (irh .ge. 40) then
          adx = adx+2.d0*ady*((-y1)+y15)
          ady = 0.
        else
          if (aw .ge. awc) then
            adawc = adawc-ady*((y40-yc)*(0.4d0-aw)/((0.4d0-awc)*(0.4d0-
     $awc)))
            ady40 = ady40+ady*(1-(0.4d0-aw)/(0.4d0-awc))
            adyc = adyc+ady*((0.4d0-aw)/(0.4d0-awc))
            ady = 0.
            adx = adx+2*adyc*y1540
            adyc = 0.
            adx = adx+2.d0*ady40*((-y140)+y1540)
            ady40 = 0.
          endif
          adx = adx+0.8d0*adawc
          adawc = 0.
        endif
      else if (x .lt. 2.d0) then
        if (irh .ge. 40) then
          adx = adx+2.d0*ady*((-y15)+y2)
          ady = 0.
        endif
      endif
      if (tso4 .gt. 0.d0) then
        adtnh4 = adtnh4+adx/tso4
        adtso4 = adtso4-adx*(tnh4/(tso4*tso4))
        adx = 0.
      endif
      admno3 = admno3+adtno3*(0.5+sign(0.5d0,mno3-0.d0))
      adtno3 = 0.
      admnh4 = admnh4+adtnh4*(0.5+sign(0.5d0,mnh4-0.d0))
      adtnh4 = 0.
      admso4 = admso4+adtso4*(0.5+sign(0.5d0,mso4-0.d0))
      adtso4 = 0.

      end SUBROUTINE ADAWATER
!-----------------------------------------------------------------------------

      subroutine adrpmares_11( in, par, adin, adout, 
     &                         I,  J,   L )

      ! References to f90 modules
      USE CHECKPT_MOD
      USE RPMARES_MOD,	ONLY : CUBIC, AWATER, ACTCOF 

#     include "CMN_SIZE"  ! Size params

C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================

C==============================================
C define parameters                            
C==============================================
      TYPE (XPLEX) floor
      parameter ( floor = xplex(1.d-30,0d0) )
      TYPE (XPLEX) mwhno3
      parameter ( mwhno3 = xplex(63.01287d0,0d0) )
      TYPE (XPLEX) mwnh3
      parameter ( mwnh3 = xplex(17.03061d0,0d0) )
      TYPE (XPLEX) mwnh4
      parameter ( mwnh4 = xplex(18.03858d0,0d0) )
      TYPE (XPLEX) mwno3
      parameter ( mwno3 = xplex(62.0049d0,0d0) )
      TYPE (XPLEX) mwso4
      parameter ( mwso4 = xplex(96.0576d0,0d0) )

C==============================================
C define common blocks
C==============================================

C==============================================
C define arguments
C==============================================
      TYPE (XPLEX) adin(5)
      TYPE (XPLEX) adout(8)
      TYPE (XPLEX) in(5)
      TYPE (XPLEX) par(2)
      INTEGER :: I, J, L

C==============================================
C define local variables
C==============================================
      TYPE (XPLEX) a0
      TYPE (XPLEX) a1
      TYPE (XPLEX) a2
      TYPE (XPLEX) ada0
      TYPE (XPLEX) ada1
      TYPE (XPLEX) ada2
      TYPE (XPLEX) adah2o
      TYPE (XPLEX) adahso4
      TYPE (XPLEX) adan(3)
      TYPE (XPLEX) adanh4
      TYPE (XPLEX) adano3
      TYPE (XPLEX) adano3_in
      TYPE (XPLEX) adaso4
      TYPE (XPLEX) adcat(2)
      TYPE (XPLEX) adcrutes(3)
      TYPE (XPLEX) aderor
      TYPE (XPLEX) aderorh
      TYPE (XPLEX) adgamahat
      TYPE (XPLEX) adgamana
      TYPE (XPLEX) adgamas1
      TYPE (XPLEX) adgamas2
      TYPE (XPLEX) adgamold
      TYPE (XPLEX) adgams(2,3)
      TYPE (XPLEX) adgnh3
      TYPE (XPLEX) adgno3
      TYPE (XPLEX) adhplus
      TYPE (XPLEX) admhso4
      TYPE (XPLEX) admna
      TYPE (XPLEX) admnh4
      TYPE (XPLEX) admso4
      TYPE (XPLEX) adrk2sa
      TYPE (XPLEX) adrkna
      TYPE (XPLEX) adrknwet
      TYPE (XPLEX) adso4
      TYPE (XPLEX) adt21
      TYPE (XPLEX) adtmasshno3
      TYPE (XPLEX) adtnh4
      TYPE (XPLEX) adtno3
      TYPE (XPLEX) adtso4
      TYPE (XPLEX) adwh2o
      TYPE (XPLEX) adxno3
      TYPE (XPLEX) adynh4
      TYPE (XPLEX) adzso4
      TYPE (XPLEX) ah2o
      TYPE (XPLEX) an(3)
      TYPE (XPLEX) anh4
      TYPE (XPLEX) ano3
      TYPE (XPLEX) cat(2)
      TYPE (XPLEX) crutes(3)
      TYPE (XPLEX) eror
      TYPE (XPLEX) erorh
      integer exit
      TYPE (XPLEX) gamaab
      TYPE (XPLEX) gamahat
      TYPE (XPLEX) gamana
      TYPE (XPLEX) gamas1
      TYPE (XPLEX) gamas2
      TYPE (XPLEX) gamas2h
      TYPE (XPLEX) gamold
      TYPE (XPLEX) gams(2,3)
      TYPE (XPLEX) gnh3
      TYPE (XPLEX) gno3
      TYPE (XPLEX) hplus
      integer ip1
      integer ip2
      integer irh
      TYPE (XPLEX) k2sa
      TYPE (XPLEX) kna
      TYPE (XPLEX) mhso4
      TYPE (XPLEX) mna
      TYPE (XPLEX) mnh4
      TYPE (XPLEX) molnu
      TYPE (XPLEX) mso4
      integer nnn
      integer nnn1
      integer nr
      TYPE (XPLEX) phibar
      TYPE (XPLEX) rh
      TYPE (XPLEX) rk2sa
      TYPE (XPLEX) rkna
      TYPE (XPLEX) rknwet
      TYPE (XPLEX) so4
      TYPE (XPLEX) t1
      TYPE (XPLEX) t2
      TYPE (XPLEX) t21
      TYPE (XPLEX) t3
      TYPE (XPLEX) t4
      TYPE (XPLEX) t6
      TYPE (XPLEX) temp
      TYPE (XPLEX) tmasshno3
      TYPE (XPLEX) tnh4
      TYPE (XPLEX) tno3
      TYPE (XPLEX) toler2
      TYPE (XPLEX) tso4
      TYPE (XPLEX) wh2o
      TYPE (XPLEX) xno3
      TYPE (XPLEX) ynh4
      TYPE (XPLEX) zso4

C----------------------------------------------
C SAVE ARGUMENTS
C----------------------------------------------
      erorh = eror

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      ada0 = 0.
      ada1 = 0.
      ada2 = 0.
      adah2o = 0.
      adahso4 = 0.
      do ip1 = 1, 3
        adan(ip1) = 0.
      end do
      adanh4 = 0.
      adano3 = 0.
      adano3_in = 0.
      adaso4 = 0.
      do ip1 = 1, 2
        adcat(ip1) = 0.
      end do
      do ip1 = 1, 3
        adcrutes(ip1) = 0.
      end do
      aderor = 0.
      adgamahat = 0.
      adgamana = 0.
      adgamas1 = 0.
      adgamas2 = 0.
      adgamold = 0.
      do ip2 = 1, 3
        do ip1 = 1, 2
          adgams(ip1,ip2) = 0.
        end do
      end do
      adgnh3 = 0.
      adgno3 = 0.
      adhplus = 0.
      admhso4 = 0.
      admna = 0.
      admnh4 = 0.
      admso4 = 0.
      adrk2sa = 0.
      adrkna = 0.
      adrknwet = 0.
      adso4 = 0.
      adt21 = 0.
      adtmasshno3 = 0.
      adtnh4 = 0.
      adtno3 = 0.
      adtso4 = 0.
      adwh2o = 0.
      adxno3 = 0.
      adynh4 = 0.
      adzso4 = 0.

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      so4 = in(1)
      gno3 = in(2)
      gnh3 = in(3)
      ano3 = in(4)
      anh4 = in(5)
      rh = par(1)
      temp = par(2)
      tso4 = max(floor,so4/mwso4)
      tno3 = max(0.d0,ano3/mwno3+gno3/mwhno3)
      tnh4 = max(0.d0,gnh3/mwnh3+anh4/mwnh4)
      tmasshno3 = max(0.d0,gno3+ano3)
      irh = nint(100.*rh)
      irh = max(1,irh)
      irh = min(99,irh)
      t6 = 8.2d-11*temp
      t1 = 298.d0/temp
      t2 = log(t1)
      t3 = t1-1.d0
      t4 = 1.d0+t2-t1
      kna = 2511000.d0*exp(29.17d0*t3+16.83d0*t4)*t6
      k2sa = 0.01015d0*exp(8.85d0*t3+25.14d0*t4)
      call awater( irh,tso4,tnh4,tno3,ah2o )
      wh2o = 0.001d0*ah2o
      zso4 = tso4/wh2o
      gamaab = 1.d0
      mnh4 = tnh4/wh2o
      ynh4 = tnh4
      adgnh3 = adgnh3+adout(8)
      adout(8) = 0.d0
      adgno3 = adgno3+adout(7)
      adout(7) = 0.d0
      adso4 = adso4+adout(6)
      adout(6) = 0.d0
      adanh4 = adanh4+adout(5)
      adout(5) = 0.d0
      adah2o = adah2o+adout(4)
      adout(4) = 0.d0
      adano3 = adano3+adout(3)
      adout(3) = 0.d0
      adahso4 = adahso4+adout(2)
      adout(2) = 0.d0
      adaso4 = adaso4+adout(1)
      adout(1) = 0.d0
      do nnn = nitr_max(I,J,L), 1, -1
        eror = erorh
        exit = 0
        call awater( irh,tso4,tnh4,tno3,ah2o )
        wh2o = 0.001d0*ah2o
        gamana = 1.d0
        gamas1 = 1.d0
        gamas2 = 1.d0
        gamold = 1.d0

        !=====================================================================
        ! CHECKPOINT
        ! The adjoint calculation needs the variables error,exit,gamana,gamas1,
        ! gamas2,gamold and wh2o at iteration nnn-1.  Rather than recompute, 
        ! use the values saved durring the forward run, but only if nnn-1 > 0
        !=====================================================================
        IF (nnn-1 .gt. 0) THEN
          eror   = eror_fwd(I,J,L,nnn-1)
          exit   = exit_fwd(I,J,L,nnn-1)
          gamana = gamana_fwd(I,J,L,nnn-1)
          gamas1 = gamas1_fwd(I,J,L,nnn-1)
          gamas2 = gamas2_fwd(I,J,L,nnn-1)
          gamold = gamold_fwd(I,J,L,nnn-1)
          wh2o   = wh2o_fwd(I,J,L,nnn-1)
        ENDIF  
  
!        do nnn1 = 1, nnn-1
!          if (exit .eq. 0) then
!            rk2sa = k2sa*gamas2*gamas2/(gamas1*gamas1*gamas1)
!            rkna = kna/(gamana*gamana)
!            rknwet = rkna*wh2o
!            t21 = zso4-mnh4
!            a2 = rk2sa+rknwet-t21
!            a1 = rk2sa*rknwet-t21*(rk2sa+rknwet)-rk2sa*zso4-rkna*tno3
!            a0 = -(t21*rk2sa*rknwet+rk2sa*rknwet*zso4+rk2sa*rkna*tno3)
!            call cubic( a2,a1,a0,nr,crutes )
!            hplus = crutes(1)
!            mso4 = rk2sa*zso4/(hplus+rk2sa)
!            mhso4 = max(1.d-10,zso4-mso4)
!            mna = rkna*tno3/(hplus+rknwet)
!            mna = max(0.,mna)
!            mna = min(mna,tno3/wh2o)
!            xno3 = mna*wh2o
!            call awater( irh,tso4,ynh4,xno3,ah2o )
!            wh2o = 0.001d0*ah2o
!            cat(1) = hplus
!            cat(2) = mnh4
!            an(1) = mso4
!            an(2) = mna
!            an(3) = mhso4
!            call actcof( cat,an,gams,molnu,phibar )
!            gamana = gams(1,2)
!            gamas1 = gams(1,1)
!            gamas2 = gams(1,3)
!            gamahat = gamas2*gamas2/(gamaab*gamaab)
!            eror = abs(gamold-gamahat)/gamold
!            gamold = gamahat
!          endif
!          if (eror .le. toler2) then
!            exit = 11
!          endif
!        end do
        gamas2h = gamas2
        if (exit .eq. 0) then
          rk2sa = k2sa*gamas2*gamas2/(gamas1*gamas1*gamas1)
          rkna = kna/(gamana*gamana)
          rknwet = rkna*wh2o
          t21 = zso4-mnh4
          a2 = rk2sa+rknwet-t21
          a1 = rk2sa*rknwet-t21*(rk2sa+rknwet)-rk2sa*zso4-rkna*tno3
          a0 = -(t21*rk2sa*rknwet+rk2sa*rknwet*zso4+rk2sa*rkna*tno3)
          call cubic( a2,a1,a0,nr,crutes )
          hplus = crutes(1)
          mso4 = rk2sa*zso4/(hplus+rk2sa)
          mhso4 = max(1.d-10,zso4-mso4)
          mna = rkna*tno3/(hplus+rknwet)
          mna = max(0.,mna)
          mna = min(mna,tno3/wh2o)
          xno3 = mna*wh2o
          ano3 = mna*wh2o*mwno3
          cat(1) = hplus
          cat(2) = mnh4
          an(1) = mso4
          an(2) = mna
          an(3) = mhso4
          call actcof( cat,an,gams,molnu,phibar )
          gamas2 = gams(1,3)
          gamahat = gamas2*gamas2/(gamaab*gamaab)
          adgamahat = adgamahat+adgamold
          adgamold = 0.
          aderorh = aderor/gamold
          adgamold = adgamold-aderor*(abs(gamold-gamahat)/(gamold*
     $gamold))
          adgamahat = adgamahat-aderorh*sign(1.d0,gamold-gamahat)
          adgamold = adgamold+aderorh*sign(1.d0,gamold-gamahat)
          aderor = 0.
          adgamas2 = adgamas2+adgamahat*(2*gamas2/(gamaab*gamaab))
          adgamahat = 0.
          adgams(1,3) = adgams(1,3)+adgamas2
          adgamas2 = 0.
          adgams(1,1) = adgams(1,1)+adgamas1
          adgamas1 = 0.
          adgams(1,2) = adgams(1,2)+adgamana
          adgamana = 0.
          call adactcof( cat,an,adcat,adan,adgams )
          admhso4 = admhso4+adan(3)
          adan(3) = 0.
          admna = admna+adan(2)
          adan(2) = 0.
          admso4 = admso4+adan(1)
          adan(1) = 0.
          admnh4 = admnh4+adcat(2)
          adcat(2) = 0.
          adhplus = adhplus+adcat(1)
          adcat(1) = 0.
          adah2o = adah2o+0.001d0*adwh2o
          adwh2o = 0.
          call adawater( irh,tso4,ynh4,xno3,adtso4,adynh4,adxno3,adah2o 
     $)
          admhso4 = admhso4+adahso4*wh2o*mwso4
          adwh2o = adwh2o+adahso4*mhso4*mwso4
          adahso4 = 0.
          admso4 = admso4+adaso4*wh2o*mwso4
          adwh2o = adwh2o+adaso4*mso4*mwso4
          adaso4 = 0.
          adano3 =adano3-adgno3*(0.5-sign(0.5d0,floor-(tmasshno3-ano3)))
          adtmasshno3 = adtmasshno3+adgno3*(0.5-sign(0.5d0,floor-
     $(tmasshno3-ano3)))
          adgno3 = 0.
          admna = admna+adano3*wh2o*mwno3
          adwh2o = adwh2o+adano3*mna*mwno3
          adano3 = 0.
          admna = admna+adxno3*wh2o
          adwh2o = adwh2o+adxno3*mna
          adxno3 = 0.
          mna = rkna*tno3/(hplus+rknwet)
          mna = max(0.,mna)
          adtno3 = adtno3+admna*((0.5-sign(0.5d0,tno3/wh2o-mna))/wh2o)
          adwh2o = adwh2o-admna*(0.5-sign(0.5d0,tno3/wh2o-mna))*(tno3/
     $(wh2o*wh2o))
          admna = admna*(0.5+sign(0.5d0,tno3/wh2o-mna))
          mna = rkna*tno3/(hplus+rknwet)
          admna = admna*(0.5-sign(0.5d0,0.-mna))
          adhplus = adhplus-admna*(rkna*tno3/((hplus+rknwet)*(hplus+
     $rknwet)))
          adrkna = adrkna+admna*(tno3/(hplus+rknwet))
          adrknwet = adrknwet-admna*(rkna*tno3/((hplus+rknwet)*(hplus+
     $rknwet)))
          adtno3 = adtno3+admna*(rkna/(hplus+rknwet))
          admna = 0.
          admso4 = admso4-admhso4*(0.5-sign(0.5d0,1.d-10-(zso4-mso4)))
          adzso4 = adzso4+admhso4*(0.5-sign(0.5d0,1.d-10-(zso4-mso4)))
          admhso4 = 0.
          adhplus = adhplus-admso4*(rk2sa*zso4/((hplus+rk2sa)*(hplus+
     $rk2sa)))
          adrk2sa = adrk2sa+admso4*(zso4/(hplus+rk2sa)-rk2sa*zso4/
     $((hplus+rk2sa)*(hplus+rk2sa)))
          adzso4 = adzso4+admso4*(rk2sa/(hplus+rk2sa))
          admso4 = 0.
          adcrutes(1) = adcrutes(1)+adhplus
          adhplus = 0.
          call adcubic( a2,a1,a0,ada2,ada1,ada0,adcrutes )
          adrk2sa = adrk2sa-ada0*(rknwet*(t21+zso4)+rkna*tno3)
          adrkna = adrkna-ada0*rk2sa*tno3
          adrknwet = adrknwet-ada0*rk2sa*(t21+zso4)
          adt21 = adt21-ada0*rk2sa*rknwet
          adtno3 = adtno3-ada0*rk2sa*rkna
          adzso4 = adzso4-ada0*rk2sa*rknwet
          ada0 = 0.
          adrk2sa = adrk2sa+ada1*(rknwet-t21-zso4)
          adrkna = adrkna-ada1*tno3
          adrknwet = adrknwet+ada1*(rk2sa-t21)
          adt21 = adt21-ada1*(rk2sa+rknwet)
          adtno3 = adtno3-ada1*rkna
          adzso4 = adzso4-ada1*rk2sa
          ada1 = 0.
          adrk2sa = adrk2sa+ada2
          adrknwet = adrknwet+ada2
          adt21 = adt21-ada2
          ada2 = 0.
          admnh4 = admnh4-adt21
          adzso4 = adzso4+adt21
          adt21 = 0.
          adrkna = adrkna+adrknwet*wh2o
          adwh2o = adwh2o+adrknwet*rkna
          adrknwet = 0.
          adgamana = adgamana-adrkna*(2*kna*gamana/(gamana*gamana*
     $gamana*gamana))
          adrkna = 0.
          gamas2 = gamas2h
          adgamas1 = adgamas1-adrk2sa*(3*k2sa*gamas2*gamas2*gamas1*
     $gamas1/(gamas1*gamas1*gamas1*gamas1*gamas1*gamas1))
          adgamas2 = adgamas2+adrk2sa*(2*k2sa*gamas2/(gamas1*gamas1*
     $gamas1))
          adrk2sa = 0.
        endif
      end do
      adtnh4 = adtnh4+adynh4
      adynh4 = 0.
      call awater( irh,tso4,tnh4,tno3,ah2o )
      wh2o = 0.001d0*ah2o
      adtnh4 = adtnh4+admnh4/wh2o
      adwh2o = adwh2o-admnh4*(tnh4/(wh2o*wh2o))
      admnh4 = 0.
      adtso4 = adtso4+adzso4/wh2o
      adwh2o = adwh2o-adzso4*(tso4/(wh2o*wh2o))
      adzso4 = 0.
      adgnh3 = 0.
      adano3 = adano3-adgno3
      adtmasshno3 = adtmasshno3+adgno3
      adgno3 = 0.
      adano3_in = adano3_in+adano3
      adano3 = 0.
      adtnh4 = adtnh4+adanh4*mwnh4
      adanh4 = 0.
      adtso4 = adtso4+adahso4*mwso4
      adahso4 = 0.
      adaso4 = 0.
      adah2o = adah2o+0.001d0*adwh2o
      adwh2o = 0.
      call adawater( irh,tso4,tnh4,tno3,adtso4,adtnh4,adtno3,adah2o )
      ano3 = in(4)
      adano3 = adano3+adtmasshno3*(0.5-sign(0.5d0,0.d0-(gno3+ano3)))
      adgno3 = adgno3+adtmasshno3*(0.5-sign(0.5d0,0.d0-(gno3+ano3)))
      adtmasshno3 = 0.
      adano3 = adano3+adano3_in
      adano3_in = 0.
      adanh4 = adanh4+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh4)
      adgnh3 = adgnh3+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh3)
      adtnh4 = 0.
      adano3 = adano3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwno3)
      adgno3 = adgno3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwhno3)
      adtno3 = 0.
      adso4 = adso4+adtso4*((0.5-sign(0.5d0,floor-so4/mwso4))/mwso4)
      adtso4 = 0.
      adin(5) = adin(5)+adanh4
      adanh4 = 0.
      adin(4) = adin(4)+adano3
      adano3 = 0.
      adin(3) = adin(3)+adgnh3
      adgnh3 = 0.
      adin(2) = adin(2)+adgno3
      adgno3 = 0.
      adin(1) = adin(1)+adso4
      adso4 = 0.

      end SUBROUTINE ADRPMARES_11
!-----------------------------------------------------------------------------


C                           DISCLAIMER
C
C   This file was generated by TAMC version 5.3.2
C
C   THE AUTHOR DOES NOT MAKE  ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR  RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
C   OR USEFULNESS  OF ANY INFORMATION OR PROCESS  DISCLOSED, OR REPRESENTS
C   THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
C   THIS CODE IS FOR NON-PROFIT-ORIENTATED ACADEMIC RESEARCH AND EDUCATION
C   ONLY.  ANY COMMERCIAL OR OTHER PROFIT-ORIENTATED USE OR  EVALUATION IS
C   STRICTLY  FORBIDDEN.  PASSING  THIS CODE  TO  ANY THIRD  PARTY IS  NOT
C   ALLOWED.
C
C   FOR COMMERCIAL OR  OTHER PROFIT-ORIENTATED APPLICATIONS PLEASE CONTACT
C   info@FastOpt.com
C
      subroutine adrpmares_12( in, par, adin, adout,
     &                         I,  J,   L )

C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================

C==============================================
C define parameters                            
C==============================================
      TYPE (XPLEX) floor
      parameter ( floor = xplex(1.d-30,0d0) )
      TYPE (XPLEX) mwhno3
      parameter ( mwhno3 = xplex(63.01287d0,0d0) )
      TYPE (XPLEX) mwnh3
      parameter ( mwnh3 = xplex(17.03061d0,0d0) )
      TYPE (XPLEX) mwnh4
      parameter ( mwnh4 = xplex(18.03858d0,0d0) )
      TYPE (XPLEX) mwno3
      parameter ( mwno3 = xplex(62.0049d0,0d0) )
      TYPE (XPLEX) mwso4
      parameter ( mwso4 = xplex(96.0576d0,0d0) )

C==============================================
C define common blocks
C==============================================
C==============================================
C define arguments
C==============================================
      TYPE (XPLEX) adin(5)
      TYPE (XPLEX) adout(8)
      TYPE (XPLEX) in(5)
      TYPE (XPLEX) par(2)
      INTEGER :: I, J, L

C==============================================
C define local variables
C==============================================
      TYPE (XPLEX) adah2o
      TYPE (XPLEX) adanh4
      TYPE (XPLEX) adano3
      TYPE (XPLEX) adano3_in
      TYPE (XPLEX) adaso4
      TYPE (XPLEX) adgnh3
      TYPE (XPLEX) adgno3
      TYPE (XPLEX) adgno3_in
      TYPE (XPLEX) adso4
      TYPE (XPLEX) adtnh4
      TYPE (XPLEX) adtno3
      TYPE (XPLEX) adtso4
      TYPE (XPLEX) anh4
      TYPE (XPLEX) ano3
      TYPE (XPLEX) gnh3
      TYPE (XPLEX) gno3
      integer irh
      TYPE (XPLEX) rh
      TYPE (XPLEX) so4
      TYPE (XPLEX) tnh4
      TYPE (XPLEX) tno3
      TYPE (XPLEX) tso4

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adah2o = 0.
      adanh4 = 0.
      adano3 = 0.
      adano3_in = 0.
      adaso4 = 0.
      adgnh3 = 0.
      adgno3 = 0.
      adgno3_in = 0.
      adso4 = 0.
      adtnh4 = 0.
      adtno3 = 0.
      adtso4 = 0.

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      so4 = in(1)
      gno3 = in(2)
      gnh3 = in(3)
      ano3 = in(4)
      anh4 = in(5)
      rh = par(1)
      tso4 = max(floor,so4/mwso4)
      tno3 = max(0.d0,ano3/mwno3+gno3/mwhno3)
      tnh4 = max(0.d0,gnh3/mwnh3+anh4/mwnh4)
      irh = nint(100.*rh)
      irh = max(1,irh)
      irh = min(99,irh)
      adgnh3 = adgnh3+adout(8)
      adout(8) = 0.d0
      adgno3 = adgno3+adout(7)
      adout(7) = 0.d0
      adso4 = adso4+adout(6)
      adout(6) = 0.d0
      adanh4 = adanh4+adout(5)
      adout(5) = 0.d0
      adah2o = adah2o+adout(4)
      adout(4) = 0.d0
      adano3 = adano3+adout(3)
      adout(3) = 0.d0
      adout(2) = 0.d0
      adaso4 = adaso4+adout(1)
      adout(1) = 0.d0
      call adawater( irh,tso4,tnh4,tno3,adtso4,adtnh4,adtno3,adah2o )
      adtso4 = adtso4+adaso4*mwso4
      adaso4 = 0.
      adano3_in = adano3_in+adano3
      adano3 = 0.
      adgno3_in = adgno3_in+adgno3
      adgno3 = 0.
      adgnh3 = 0.
      adtnh4 = adtnh4+adanh4*mwnh4
      adanh4 = 0.
      adano3 = adano3+adano3_in
      adano3_in = 0.
      adgno3 = adgno3+adgno3_in
      adgno3_in = 0.
      adanh4 = adanh4+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh4)
      adgnh3 = adgnh3+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh3)
      adtnh4 = 0.
      adano3 = adano3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwno3)
      adgno3 = adgno3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwhno3)
      adtno3 = 0.
      adso4 = adso4+adtso4*((0.5-sign(0.5d0,floor-so4/mwso4))/mwso4)
      adtso4 = 0.
      adin(5) = adin(5)+adanh4
      adanh4 = 0.
      adin(4) = adin(4)+adano3
      adano3 = 0.
      adin(3) = adin(3)+adgnh3
      adgnh3 = 0.
      adin(2) = adin(2)+adgno3
      adgno3 = 0.
      adin(1) = adin(1)+adso4
      adso4 = 0.

      end SUBROUTINE ADRPMARES_12


C                           DISCLAIMER
C
C   This file was generated by TAMC version 5.3.2
C
C   THE AUTHOR DOES NOT MAKE  ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR  RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
C   OR USEFULNESS  OF ANY INFORMATION OR PROCESS  DISCLOSED, OR REPRESENTS
C   THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
C   THIS CODE IS FOR NON-PROFIT-ORIENTATED ACADEMIC RESEARCH AND EDUCATION
C   ONLY.  ANY COMMERCIAL OR OTHER PROFIT-ORIENTATED USE OR  EVALUATION IS
C   STRICTLY  FORBIDDEN.  PASSING  THIS CODE  TO  ANY THIRD  PARTY IS  NOT
C   ALLOWED.
C
C   FOR COMMERCIAL OR  OTHER PROFIT-ORIENTATED APPLICATIONS PLEASE CONTACT
C   info@FastOpt.com
C
      subroutine adrpmares_2( in, par, adin, adout,
     &                         I,  J,   L )

C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================

C==============================================
C define parameters                            
C==============================================
      TYPE (XPLEX) floor
      parameter ( floor = xplex(1.d-30,0d0) )
      TYPE (XPLEX) mwno3
      parameter ( mwno3 = xplex(62.0049d0,0d0) )
      TYPE (XPLEX) minno3
      parameter ( minno3 = xplex(1.d-6/mwno3%r,0d0) )
      TYPE (XPLEX) mwso4
      parameter ( mwso4 = xplex(96.0576d0,0d0) )
      TYPE (XPLEX) minso4
      parameter ( minso4 = xplex(1.d-6/mwso4%r,0d0) )
      TYPE (XPLEX) mwhno3
      parameter ( mwhno3 = xplex(63.01287d0,0d0) )

C==============================================
C define arguments
C==============================================
      TYPE (XPLEX) adin(5)
      TYPE (XPLEX) adout(8)
      TYPE (XPLEX) in(5)
      TYPE (XPLEX) par(2)
      INTEGER :: I, J, L

C==============================================
C define local variables
C==============================================
      TYPE (XPLEX) adanh4
      TYPE (XPLEX) adano3
      TYPE (XPLEX) adaso4
      TYPE (XPLEX) adgnh3
      TYPE (XPLEX) adgno3
      TYPE (XPLEX) adso4
      TYPE (XPLEX) adtno3
      TYPE (XPLEX) adtso4
      TYPE (XPLEX) anh4
      TYPE (XPLEX) ano3
      TYPE (XPLEX) aso4
      integer exit
      TYPE (XPLEX) gnh3
      TYPE (XPLEX) gno3
      TYPE (XPLEX) rh
      TYPE (XPLEX) so4
      TYPE (XPLEX) tno3
      TYPE (XPLEX) tso4

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adanh4 = 0.
      adano3 = 0.
      adaso4 = 0.
      adgnh3 = 0.
      adgno3 = 0.
      adso4 = 0.
      adtno3 = 0.
      adtso4 = 0.

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      so4 = in(1)
      gno3 = in(2)
      gnh3 = in(3)
      ano3 = in(4)
      anh4 = in(5)
      rh = par(1)
      aso4 = 0.d0
      exit = 0
      if (rh .lt. 0.01) then
        exit = 1
      endif
      if (exit .eq. 0) then
        tso4 = max(floor,so4/mwso4)
        aso4 = so4
        tno3 = max(0.d0,ano3/mwno3+gno3/mwhno3)
      endif
      adgnh3 = adgnh3+adout(8)
      adout(8) = 0.d0
      adgno3 = adgno3+adout(7)
      adout(7) = 0.d0
      adso4 = adso4+adout(6)
      adout(6) = 0.d0
      adanh4 = adanh4+adout(5)
      adout(5) = 0.d0
      adout(4) = 0.d0
      adano3 = adano3+adout(3)
      adout(3) = 0.d0
      adout(2) = 0.d0
      adaso4 = adaso4+adout(1)
      adout(1) = 0.d0
      if (tso4 .lt. minso4 .and. tno3 .lt. minno3) then
        if (exit .eq. 0) then
          adgno3 = adgno3*(0.5-sign(0.5d0,floor-gno3))
          adgnh3 = adgnh3*(0.5-sign(0.5d0,floor-gnh3))
          adanh4 = adanh4*(0.5-sign(0.5d0,floor-anh4))
          adano3 = adano3*(0.5-sign(0.5d0,floor-ano3))
          adaso4 = adaso4*(0.5-sign(0.5d0,floor-aso4))
        endif
      endif
      if (exit .eq. 0) then
        adano3 = adano3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwno3)
        adgno3 = adgno3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwhno3)
        adtno3 = 0.
        adso4 = adso4+adaso4
        adaso4 = 0.
        adso4 = adso4+adtso4*((0.5-sign(0.5d0,floor-so4/mwso4))/mwso4)
        adtso4 = 0.
      endif
      adin(5) = adin(5)+adanh4
      adanh4 = 0.
      adin(4) = adin(4)+adano3
      adano3 = 0.
      adin(3) = adin(3)+adgnh3
      adgnh3 = 0.
      adin(2) = adin(2)+adgno3
      adgno3 = 0.
      adin(1) = adin(1)+adso4
      adso4 = 0.

      end SUBROUTINE ADRPMARES_2
!-----------------------------------------------------------------------------

C                           DISCLAIMER
C
C   This file was generated by TAMC version 5.3.2
C
C   THE AUTHOR DOES NOT MAKE  ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR  RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
C   OR USEFULNESS  OF ANY INFORMATION OR PROCESS  DISCLOSED, OR REPRESENTS
C   THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
C   THIS CODE IS FOR NON-PROFIT-ORIENTATED ACADEMIC RESEARCH AND EDUCATION
C   ONLY.  ANY COMMERCIAL OR OTHER PROFIT-ORIENTATED USE OR  EVALUATION IS
C   STRICTLY  FORBIDDEN.  PASSING  THIS CODE  TO  ANY THIRD  PARTY IS  NOT
C   ALLOWED.
C
C   FOR COMMERCIAL OR  OTHER PROFIT-ORIENTATED APPLICATIONS PLEASE CONTACT
C   info@FastOpt.com
C
      subroutine adrpmares_3( in, par, adin, adout,
     &                         I,  J,   L )

C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================

C==============================================
C define parameters                            
C==============================================
      TYPE (XPLEX) floor
      parameter ( floor = xplex(1.d-30,0d0) )
      TYPE (XPLEX) mwhno3
      parameter ( mwhno3 = xplex(63.01287d0,0d0) )
      TYPE (XPLEX) mwnh3
      parameter ( mwnh3 = xplex(17.03061d0,0d0) )
      TYPE (XPLEX) mwnh4
      parameter ( mwnh4 = xplex(18.03858d0,0d0) )
      TYPE (XPLEX) mwno3
      parameter ( mwno3 = xplex(62.0049d0,0d0) )
      TYPE (XPLEX) mwso4
      parameter ( mwso4 = xplex(96.0576d0,0d0) )

C==============================================
C define common blocks
C==============================================
C==============================================
C define arguments
C==============================================
      TYPE (XPLEX) adin(5)
      TYPE (XPLEX) adout(8)
      TYPE (XPLEX) in(5)
      TYPE (XPLEX) par(2)
      INTEGER :: I, J, L 

C==============================================
C define local variables
C==============================================
      TYPE (XPLEX) adah2o
      TYPE (XPLEX) adanh4
      TYPE (XPLEX) adano3
      TYPE (XPLEX) adano3_in
      TYPE (XPLEX) adaso4
      TYPE (XPLEX) adgnh3
      TYPE (XPLEX) adgno3
      TYPE (XPLEX) adgno3_in
      TYPE (XPLEX) adso4
      TYPE (XPLEX) adtnh4
      TYPE (XPLEX) adtno3
      TYPE (XPLEX) adtso4
      TYPE (XPLEX) adtwoso4
      TYPE (XPLEX) adwh2o
      TYPE (XPLEX) adynh4
      TYPE (XPLEX) anh4
      TYPE (XPLEX) ano3
      TYPE (XPLEX) gnh3
      TYPE (XPLEX) gno3
      integer irh
      TYPE (XPLEX) rh
      TYPE (XPLEX) so4
      TYPE (XPLEX) tnh4
      TYPE (XPLEX) tno3
      TYPE (XPLEX) tso4
      TYPE (XPLEX) twoso4
      TYPE (XPLEX) ynh4

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adah2o = 0.
      adanh4 = 0.
      adano3 = 0.
      adano3_in = 0.
      adaso4 = 0.
      adgnh3 = 0.
      adgno3 = 0.
      adgno3_in = 0.
      adso4 = 0.
      adtnh4 = 0.
      adtno3 = 0.
      adtso4 = 0.
      adtwoso4 = 0.
      adwh2o = 0.
      adynh4 = 0.

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      so4 = in(1)
      gno3 = in(2)
      gnh3 = in(3)
      ano3 = in(4)
      anh4 = in(5)
      rh = par(1)
      tso4 = max(floor,so4/mwso4)
      tno3 = max(0.d0,ano3/mwno3+gno3/mwhno3)
      tnh4 = max(0.d0,gnh3/mwnh3+anh4/mwnh4)
      irh = nint(100.*rh)
      irh = max(1,irh)
      irh = min(99,irh)
      twoso4 = 2.*tso4
      ynh4 = twoso4
      adgnh3 = adgnh3+adout(8)
      adout(8) = 0.d0
      adgno3 = adgno3+adout(7)
      adout(7) = 0.d0
      adso4 = adso4+adout(6)
      adout(6) = 0.d0
      adanh4 = adanh4+adout(5)
      adout(5) = 0.d0
      adah2o = adah2o+adout(4)
      adout(4) = 0.d0
      adano3 = adano3+adout(3)
      adout(3) = 0.d0
      adout(2) = 0.d0
      adaso4 = adaso4+adout(1)
      adout(1) = 0.d0
      adano3_in = adano3_in+adano3
      adano3 = 0.
      adgno3_in = adgno3_in+adgno3
      adgno3 = 0.
      adtnh4 = adtnh4+adgnh3*mwnh3*(0.5-sign(0.5d0,floor-(tnh4-ynh4)))
      adynh4 = adynh4-adgnh3*mwnh3*(0.5-sign(0.5d0,floor-(tnh4-ynh4)))
      adgnh3 = 0.
      adynh4 = adynh4+adanh4*mwnh4
      adanh4 = 0.
      adtso4 = adtso4+adaso4*mwso4
      adaso4 = 0.
      adtwoso4 = adtwoso4+adynh4
      adynh4 = 0.
      adwh2o = adwh2o+1000*adah2o
      adah2o = 0.
      adah2o = adah2o+0.001d0*adwh2o
      adwh2o = 0.
      ynh4 = twoso4
      call adawater( irh,tso4,ynh4,tno3,adtso4,adynh4,adtno3,adah2o )
      adtwoso4 = adtwoso4+adynh4
      adynh4 = 0.
      adtso4 = adtso4+2*adtwoso4
      adtwoso4 = 0.
      adano3 = adano3+adano3_in
      adano3_in = 0.
      adgno3 = adgno3+adgno3_in
      adgno3_in = 0.
      adanh4 = adanh4+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh4)
      adgnh3 = adgnh3+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh3)
      adtnh4 = 0.
      adano3 = adano3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwno3)
      adgno3 = adgno3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwhno3)
      adtno3 = 0.
      adso4 = adso4+adtso4*((0.5-sign(0.5d0,floor-so4/mwso4))/mwso4)
      adtso4 = 0.
      adin(5) = adin(5)+adanh4
      adanh4 = 0.
      adin(4) = adin(4)+adano3
      adano3 = 0.
      adin(3) = adin(3)+adgnh3
      adgnh3 = 0.
      adin(2) = adin(2)+adgno3
      adgno3 = 0.
      adin(1) = adin(1)+adso4
      adso4 = 0.

      end SUBROUTINE ADRPMARES_3

!-----------------------------------------------------------------------------


C                           DISCLAIMER
C
C   This file was generated by TAMC version 5.3.2
C
C   THE AUTHOR DOES NOT MAKE  ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR  RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
C   OR USEFULNESS  OF ANY INFORMATION OR PROCESS  DISCLOSED, OR REPRESENTS
C   THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
C   THIS CODE IS FOR NON-PROFIT-ORIENTATED ACADEMIC RESEARCH AND EDUCATION
C   ONLY.  ANY COMMERCIAL OR OTHER PROFIT-ORIENTATED USE OR  EVALUATION IS
C   STRICTLY  FORBIDDEN.  PASSING  THIS CODE  TO  ANY THIRD  PARTY IS  NOT
C   ALLOWED.
C
C   FOR COMMERCIAL OR  OTHER PROFIT-ORIENTATED APPLICATIONS PLEASE CONTACT
C   info@FastOpt.com
C
      subroutine adrpmares_4( in, par, adin, adout, 
     &                         I,  J,   L )

C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================

C==============================================
C define parameters                            
C==============================================
      TYPE (XPLEX) floor
      parameter ( floor = xplex(1.d-30,0d0) )
      TYPE (XPLEX) mwhno3
      parameter ( mwhno3 = xplex(63.01287d0,0d0) )
      TYPE (XPLEX) mwnh3
      parameter ( mwnh3 = xplex(17.03061d0,0d0) )
      TYPE (XPLEX) mwnh4
      parameter ( mwnh4 = xplex(18.03858d0,0d0) )
      TYPE (XPLEX) mwno3
      parameter ( mwno3 = xplex(62.0049d0,0d0) )
      TYPE (XPLEX) mwso4
      parameter ( mwso4 = xplex(96.0576d0,0d0) )

C==============================================
C define common blocks
C==============================================
C==============================================
C define arguments
C==============================================
      TYPE (XPLEX) adin(5)
      TYPE (XPLEX) adout(8)
      TYPE (XPLEX) in(5)
      TYPE (XPLEX) par(2)
      INTEGER :: I, J, L

C==============================================
C define local variables
C==============================================
      TYPE (XPLEX) aa
      TYPE (XPLEX) adah2o
      TYPE (XPLEX) adanh4
      TYPE (XPLEX) adano3
      TYPE (XPLEX) adaso4
      TYPE (XPLEX) adbb
      TYPE (XPLEX) adcc
      TYPE (XPLEX) addd
      TYPE (XPLEX) addisc
      TYPE (XPLEX) adfnh3
      TYPE (XPLEX) adgnh3
      TYPE (XPLEX) adgno3
      TYPE (XPLEX) adso4
      TYPE (XPLEX) adtmasshno3
      TYPE (XPLEX) adtnh4
      TYPE (XPLEX) adtno3
      TYPE (XPLEX) adtso4
      TYPE (XPLEX) adtwoso4
      TYPE (XPLEX) adwh2o
      TYPE (XPLEX) adxno3
      TYPE (XPLEX) adxxq
      TYPE (XPLEX) adynh4
      TYPE (XPLEX) anh4
      TYPE (XPLEX) ano3
      TYPE (XPLEX) bb
      TYPE (XPLEX) cc
      TYPE (XPLEX) convt
      TYPE (XPLEX) dd
      TYPE (XPLEX) disc
      TYPE (XPLEX) fnh3
      TYPE (XPLEX) gnh3
      TYPE (XPLEX) gno3
      integer irh
      TYPE (XPLEX) k3
      TYPE (XPLEX) rh
      TYPE (XPLEX) so4
      TYPE (XPLEX) temp
      TYPE (XPLEX) tmasshno3
      TYPE (XPLEX) tnh4
      TYPE (XPLEX) tno3
      TYPE (XPLEX) tso4
      TYPE (XPLEX) twoso4
      TYPE (XPLEX) xno3
      TYPE (XPLEX) xxq
      TYPE (XPLEX) ynh4

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adah2o = 0.
      adanh4 = 0.
      adano3 = 0.
      adaso4 = 0.
      adbb = 0.
      adcc = 0.
      addd = 0.
      addisc = 0.
      adfnh3 = 0.
      adgnh3 = 0.
      adgno3 = 0.
      adso4 = 0.
      adtmasshno3 = 0.
      adtnh4 = 0.
      adtno3 = 0.
      adtso4 = 0.
      adtwoso4 = 0.
      adwh2o = 0.
      adxno3 = 0.
      adxxq = 0.
      adynh4 = 0.

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      so4 = in(1)
      gno3 = in(2)
      gnh3 = in(3)
      ano3 = in(4)
      anh4 = in(5)
      rh = par(1)
      temp = par(2)
      tso4 = max(floor,so4/mwso4)
      tno3 = max(0.d0,ano3/mwno3+gno3/mwhno3)
      tnh4 = max(0.d0,gnh3/mwnh3+anh4/mwnh4)
      tmasshno3 = max(0.d0,gno3+ano3)
      irh = nint(100.*rh)
      irh = max(1,irh)
      irh = min(99,irh)
      convt = 1.d0/(0.082d0*temp)
      k3 = exp(118.87-24084./temp-6.025*log(temp))
      k3 = k3*convt*convt
      twoso4 = 2.*tso4
      fnh3 = tnh4-twoso4
      cc = tno3*fnh3-k3
      if (cc .le. 0.d0) then
        xno3 = 0.d0
      else
        aa = 1.d0
        bb = -(tno3+fnh3)
        disc = bb*bb-4.*cc
        dd = sqrt(disc)
        xxq = -(0.5d0*(bb+sign(1.d0,bb)*dd))
        xno3 = min(xxq/aa,cc/xxq)
      endif
      ynh4 = twoso4+xno3
      ano3 = xno3*mwno3
      adgnh3 = adgnh3+adout(8)
      adout(8) = 0.d0
      adgno3 = adgno3+adout(7)
      adout(7) = 0.d0
      adso4 = adso4+adout(6)
      adout(6) = 0.d0
      adanh4 = adanh4+adout(5)
      adout(5) = 0.d0
      adah2o = adah2o+adout(4)
      adout(4) = 0.d0
      adano3 = adano3+adout(3)
      adout(3) = 0.d0
      adout(2) = 0.d0
      adaso4 = adaso4+adout(1)
      adout(1) = 0.d0
      adano3 = adano3-adgno3*(0.5-sign(0.5d0,floor-(tmasshno3-ano3)))
      adtmasshno3 = adtmasshno3+adgno3*(0.5-sign(0.5d0,floor-(tmasshno3-
     $ano3)))
      adgno3 = 0.
      adtnh4 = adtnh4+adgnh3*mwnh3*(0.5-sign(0.5d0,floor-(tnh4-ynh4)))
      adynh4 = adynh4-adgnh3*mwnh3*(0.5-sign(0.5d0,floor-(tnh4-ynh4)))
      adgnh3 = 0.
      adynh4 = adynh4+adanh4*mwnh4
      adanh4 = 0.
      adxno3 = adxno3+adano3*mwno3
      adano3 = 0.
      adtso4 = adtso4+adaso4*mwso4
      adaso4 = 0.
      adtwoso4 = adtwoso4+adynh4
      adxno3 = adxno3+adynh4
      adynh4 = 0.
      adwh2o = adwh2o+1000*adah2o
      adah2o = 0.
      if (cc .le. 0.d0) then
      else
        adcc = adcc+adxno3*((0.5-sign(0.5d0,cc/xxq-xxq/aa))/xxq)
        adxxq = adxxq+adxno3*((0.5+sign(0.5d0,cc/xxq-xxq/aa))/aa-(0.5-
     $sign(0.5d0,cc/xxq-xxq/aa))*(cc/(xxq*xxq)))
        adxno3 = 0.
        adbb = adbb-0.5d0*adxxq
        addd = addd-0.5d0*adxxq*sign(1.d0,bb)
        adxxq = 0.
        addisc = addisc+addd*(1./(2.*sqrt(disc)))
        addd = 0.
        adbb = adbb+2*addisc*bb
        adcc = adcc-4*addisc
        addisc = 0.
        adfnh3 = adfnh3-adbb
        adtno3 = adtno3-adbb
        adbb = 0.
      endif
      adfnh3 = adfnh3+adcc*tno3
      adtno3 = adtno3+adcc*fnh3
      adcc = 0.
      adtnh4 = adtnh4+adfnh3
      adtwoso4 = adtwoso4-adfnh3
      adfnh3 = 0.
      adah2o = adah2o+0.001d0*adwh2o
      adwh2o = 0.
      ynh4 = twoso4
      call adawater( irh,tso4,ynh4,tno3,adtso4,adynh4,adtno3,adah2o )
      adtwoso4 = adtwoso4+adynh4
      adynh4 = 0.
      adtso4 = adtso4+2*adtwoso4
      adtwoso4 = 0.
      ano3 = in(4)
      adano3 = adano3+adtmasshno3*(0.5-sign(0.5d0,0.d0-(gno3+ano3)))
      adgno3 = adgno3+adtmasshno3*(0.5-sign(0.5d0,0.d0-(gno3+ano3)))
      adtmasshno3 = 0.
      adanh4 = adanh4+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh4)
      adgnh3 = adgnh3+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh3)
      adtnh4 = 0.
      adano3 = adano3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwno3)
      adgno3 = adgno3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwhno3)
      adtno3 = 0.
      adso4 = adso4+adtso4*((0.5-sign(0.5d0,floor-so4/mwso4))/mwso4)
      adtso4 = 0.
      adin(5) = adin(5)+adanh4
      adanh4 = 0.
      adin(4) = adin(4)+adano3
      adano3 = 0.
      adin(3) = adin(3)+adgnh3
      adgnh3 = 0.
      adin(2) = adin(2)+adgno3
      adgno3 = 0.
      adin(1) = adin(1)+adso4
      adso4 = 0.

      end SUBROUTINE ADRPMARES_4

!-----------------------------------------------------------------------------


C                           DISCLAIMER
C
C   This file was generated by TAMC version 5.3.2
C
C   THE AUTHOR DOES NOT MAKE  ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR  RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
C   OR USEFULNESS  OF ANY INFORMATION OR PROCESS  DISCLOSED, OR REPRESENTS
C   THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
C   THIS CODE IS FOR NON-PROFIT-ORIENTATED ACADEMIC RESEARCH AND EDUCATION
C   ONLY.  ANY COMMERCIAL OR OTHER PROFIT-ORIENTATED USE OR  EVALUATION IS
C   STRICTLY  FORBIDDEN.  PASSING  THIS CODE  TO  ANY THIRD  PARTY IS  NOT
C   ALLOWED.
C
C   FOR COMMERCIAL OR  OTHER PROFIT-ORIENTATED APPLICATIONS PLEASE CONTACT
C   info@FastOpt.com
C
      subroutine adrpmares_6( in, par, adin, adout,
     &                         I,  J,   L )

      ! References to f90 modules
      USE CHECKPT_MOD
      USE RPMARES_MOD,  	ONLY : AWATER, ACTCOF

#     include "CMN_SIZE"  ! Size params
 
C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================

C==============================================
C define parameters                            
C==============================================
      TYPE (XPLEX) floor
      parameter ( floor = xplex(1.d-30,0d0) )
      TYPE (XPLEX) mwhno3
      parameter ( mwhno3 = xplex(63.01287d0,0d0) )
      TYPE (XPLEX) mwnh3
      parameter ( mwnh3 = xplex(17.03061d0,0d0) )
      TYPE (XPLEX) mwnh4
      parameter ( mwnh4 = xplex(18.03858d0,0d0) )
      TYPE (XPLEX) mwno3
      parameter ( mwno3 = xplex(62.0049d0,0d0) )
      TYPE (XPLEX) mwso4
      parameter ( mwso4 = xplex(96.0576d0,0d0) )

C==============================================
C define common blocks
C==============================================

C==============================================
C define arguments
C============================================= =
      TYPE (XPLEX) adin(5)
      TYPE (XPLEX) adout(8)
      TYPE (XPLEX) in(5)
      TYPE (XPLEX) par(2)
      INTEGER :: I, J, L

C==============================================
C define local variables
C==============================================
      TYPE (XPLEX) aa
      TYPE (XPLEX) adaa
      TYPE (XPLEX) adah2o
      TYPE (XPLEX) adan(3)
      TYPE (XPLEX) adanh4
      TYPE (XPLEX) adano3
      TYPE (XPLEX) adaso4
      TYPE (XPLEX) adbb
      TYPE (XPLEX) adcat(2)
      TYPE (XPLEX) adcc
      TYPE (XPLEX) addd
      TYPE (XPLEX) addisc
      TYPE (XPLEX) aderor
      TYPE (XPLEX) aderorh
      TYPE (XPLEX) adgamaan
      TYPE (XPLEX) adgamold
      TYPE (XPLEX) adgams(2,3)
      TYPE (XPLEX) adgasqd
      TYPE (XPLEX) adgnh3
      TYPE (XPLEX) adgno3
      TYPE (XPLEX) adkw2
      TYPE (XPLEX) adman
      TYPE (XPLEX) admas
      TYPE (XPLEX) admnh4
      TYPE (XPLEX) adrr1
      TYPE (XPLEX) adrr2
      TYPE (XPLEX) adso4
      TYPE (XPLEX) adtmasshno3
      TYPE (XPLEX) adtnh4
      TYPE (XPLEX) adtno3
      TYPE (XPLEX) adtso4
      TYPE (XPLEX) adtwoso4
      TYPE (XPLEX) adwh2o
      TYPE (XPLEX) adwsqd
      TYPE (XPLEX) adxno3
      TYPE (XPLEX) adxxq
      TYPE (XPLEX) adynh4
      TYPE (XPLEX) ah2o
      TYPE (XPLEX) an(3)
      TYPE (XPLEX) anh4
      TYPE (XPLEX) ano3
      TYPE (XPLEX) bb
      TYPE (XPLEX) cat(2)
      TYPE (XPLEX) cc
      TYPE (XPLEX) dd
      TYPE (XPLEX) disc
      TYPE (XPLEX) eror
      integer exit
      TYPE (XPLEX) gamaan
      TYPE (XPLEX) gamaanh
      TYPE (XPLEX) gamold
      TYPE (XPLEX) gams(2,3)
      TYPE (XPLEX) gasqd
      TYPE (XPLEX) gnh3
      TYPE (XPLEX) gno3
      integer ip1
      integer ip2
      integer irh
      TYPE (XPLEX) k1a
      TYPE (XPLEX) kan
      TYPE (XPLEX) khat
      TYPE (XPLEX) kna
      TYPE (XPLEX) kph
      TYPE (XPLEX) kw
      TYPE (XPLEX) kw2
      TYPE (XPLEX) man
      TYPE (XPLEX) mas
      TYPE (XPLEX) mnh4
      TYPE (XPLEX) molnu
      integer nnn
      integer nnn1
      TYPE (XPLEX) phibar
      TYPE (XPLEX) rh
      TYPE (XPLEX) rr1
      TYPE (XPLEX) rr2
      TYPE (XPLEX) so4
      TYPE (XPLEX) t1
      TYPE (XPLEX) t2
      TYPE (XPLEX) t3
      TYPE (XPLEX) t4
      TYPE (XPLEX) t6
      TYPE (XPLEX) temp
      TYPE (XPLEX) tmasshno3
      TYPE (XPLEX) tnh4
      TYPE (XPLEX) tno3
      TYPE (XPLEX) toler1
      TYPE (XPLEX) tso4
      TYPE (XPLEX) twoso4
      TYPE (XPLEX) wh2o
      TYPE (XPLEX) wh2oh
      TYPE (XPLEX) wsqd
      TYPE (XPLEX) xno3
      TYPE (XPLEX) xxq
      TYPE (XPLEX) ynh4
      TYPE (XPLEX) ynh4h

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adaa = 0.
      adah2o = 0.
      do ip1 = 1, 3
        adan(ip1) = 0.
      end do
      adanh4 = 0.
      adano3 = 0.
      adaso4 = 0.
      adbb = 0.
      do ip1 = 1, 2
        adcat(ip1) = 0.
      end do
      adcc = 0.
      addd = 0.
      addisc = 0.
      aderor = 0.
      adgamaan = 0.
      adgamold = 0.
      do ip2 = 1, 3
        do ip1 = 1, 2
          adgams(ip1,ip2) = 0.
        end do
      end do
      adgasqd = 0.
      adgnh3 = 0.
      adgno3 = 0.
      adkw2 = 0.
      adman = 0.
      admas = 0.
      admnh4 = 0.
      adrr1 = 0.
      adrr2 = 0.
      adso4 = 0.
      adtmasshno3 = 0.
      adtnh4 = 0.
      adtno3 = 0.
      adtso4 = 0.
      adtwoso4 = 0.
      adwh2o = 0.
      adwsqd = 0.
      adxno3 = 0.
      adxxq = 0.
      adynh4 = 0.

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      so4 = in(1)
      gno3 = in(2)
      gnh3 = in(3)
      ano3 = in(4)
      anh4 = in(5)
      rh = par(1)
      temp = par(2)
      tso4 = max(floor,so4/mwso4)
      tno3 = max(0.d0,ano3/mwno3+gno3/mwhno3)
      tnh4 = max(0.d0,gnh3/mwnh3+anh4/mwnh4)
      tmasshno3 = max(0.d0,gno3+ano3)
      irh = nint(100.*rh)
      irh = max(1,irh)
      irh = min(99,irh)
      t6 = 8.2d-11*temp
      t1 = 298.d0/temp
      t2 = log(t1)
      t3 = t1-1.d0
      t4 = 1.d0+t2-t1
      kna = 2511000.d0*exp(29.17d0*t3+16.83d0*t4)*t6
      k1a = 0.00001805d0*exp((-(1.5d0*t3))+26.92d0*t4)
      kw = 1.01d-14*exp((-(22.52d0*t3))+26.92d0*t4)
      kph = 57.639d0*exp(13.79d0*t3-5.39d0*t4)*t6
      khat = kph*k1a/kw
      kan = kna*khat
      toler1 = 0.00001d0
      twoso4 = 2.*tso4
      adgnh3 = adgnh3+adout(8)
      adout(8) = 0.d0
      adgno3 = adgno3+adout(7)
      adout(7) = 0.d0
      adso4 = adso4+adout(6)
      adout(6) = 0.d0
      adanh4 = adanh4+adout(5)
      adout(5) = 0.d0
      adah2o = adah2o+adout(4)
      adout(4) = 0.d0
      adano3 = adano3+adout(3)
      adout(3) = 0.d0
      adout(2) = 0.d0
      adaso4 = adaso4+adout(1)
      adout(1) = 0.d0
      !do nnn = 50, 1, -1
      do nnn = nitr_max(I,J,L), 1, -1
        exit = 0
        gamold = 1.d0
        gamaan = 0.1
        ynh4 = twoso4
        call awater( irh,tso4,ynh4,tno3,ah2o )
        wh2o = 0.001d0*ah2o
        ynh4 = twoso4

        !=====================================================================
        ! CHECKPOINT 
        ! The do loop below was soley for the recomputation of ynh4h, wh2oh,
        ! gamaanhm, exit and gamold.
        ! Instead of recomputing these quantities, get them from the forward
        ! calculation, except if nnn is 1, in which case they are just their
        ! normal values
        !===================================================================== 
        IF (nnn-1 .gt. 0) THEN
          ynh4h = ynh4_fwd(I,J,L,nnn-1)
          wh2oh = wh2o_fwd(I,J,L,nnn-1)
          gamaanh = gamaan_fwd(I,J,L,nnn-1)
          gamold = gamold_fwd(I,J,L,nnn-1)
          exit = exit_fwd(I,J,L,nnn-1)
        ENDIF

!        do nnn1 = 1, nnn-1
!          gasqd = gamaan*gamaan
!          wsqd = wh2o*wh2o
!          kw2 = kan*wsqd/gasqd
!          aa = 1.-kw2
!          bb = twoso4+kw2*(tno3+tnh4-twoso4)
!          cc = -(kw2*tno3*(tnh4-twoso4))
!          disc = bb*bb-4.*aa*cc
!          if (aa .ne. 0.d0) then
!            dd = sqrt(disc)
!            xxq = -(0.5d0*(bb+sign(1.d0,bb)*dd))
!            rr1 = xxq/aa
!            rr2 = cc/xxq
!            if (rr1*rr2 .lt. 0.d0) then
!              xno3 = max(rr1,rr2)
!            else
!              xno3 = min(rr1,rr2)
!            endif
!          else
!            xno3 = -(cc/bb)
!          endif
!          xno3 = min(xno3,tno3)
!          call awater( irh,tso4,ynh4,xno3,ah2o )
!          wh2o = 0.001*ah2o
!          man = xno3/wh2o
!          mas = tso4/wh2o
!          mnh4 = 2.*mas+man
!          ynh4 = mnh4*wh2o
!          cat(1) = 0.
!          cat(2) = mnh4
!          an(1) = mas
!          an(2) = man
!          an(3) = 0.
!          call actcof( cat,an,gams,molnu,phibar )
!          gamaan = gams(2,2)
!          eror = abs(gamold-gamaan)/gamold
!          gamold = gamaan
!          if (eror .le. toler1) then
!            if (exit .eq. 0) then
!              exit = 6
!            endif
!          endif
!        end do
        ynh4h = ynh4
        wh2oh = wh2o
        gamaanh = gamaan
        gasqd = gamaan*gamaan
        wsqd = wh2o*wh2o
        kw2 = kan*wsqd/gasqd
        aa = 1.-kw2
        bb = twoso4+kw2*(tno3+tnh4-twoso4)
        cc = -(kw2*tno3*(tnh4-twoso4))
        disc = bb*bb-4.*aa*cc
        if (aa .ne. 0.d0) then
          dd = sqrt(disc)
          xxq = -(0.5d0*(bb+sign(1.d0,bb)*dd))
          rr1 = xxq/aa
          rr2 = cc/xxq
          if (rr1*rr2 .lt. 0.d0) then
            xno3 = max(rr1,rr2)
          else
            xno3 = min(rr1,rr2)
          endif
        else
          xno3 = -(cc/bb)
        endif
        xno3 = min(xno3,tno3)
        call awater( irh,tso4,ynh4,xno3,ah2o )
        wh2o = 0.001*ah2o
        man = xno3/wh2o
        mas = tso4/wh2o
        mnh4 = 2.*mas+man
        ynh4 = mnh4*wh2o
        cat(1) = 0.
        cat(2) = mnh4
        an(1) = mas
        an(2) = man
        an(3) = 0.
        call actcof( cat,an,gams,molnu,phibar )
        gamaan = gams(2,2)
        eror = abs(gamold-gamaan)/gamold
        if (eror .le. toler1) then
          if (exit .eq. 0) then
            ano3 = xno3*mwno3
            adwh2o = adwh2o+1000*adah2o
            adah2o = 0.
          adtnh4 = adtnh4+adgnh3*mwnh3*(0.5-sign(0.5d0,floor-(tnh4-ynh4)
     $))
          adynh4 = adynh4-adgnh3*mwnh3*(0.5-sign(0.5d0,floor-(tnh4-ynh4)
     $))
            adgnh3 = 0.
          adano3 = adano3-adgno3*(0.5-sign(0.5d0,floor-(tmasshno3-ano3))
     $)
            adtmasshno3 = adtmasshno3+adgno3*(0.5d0-sign(0.5d0,floor-
     $(tmasshno3-ano3)))
            adgno3 = 0.
            adynh4 = adynh4+adanh4*mwnh4
            adanh4 = 0.
            adxno3 = adxno3+adano3*mwno3
            adano3 = 0.
            adtso4 = adtso4+adaso4*mwso4
            adaso4 = 0.
          endif
        endif
        adgamaan = adgamaan+adgamold
        adgamold = 0.
        aderorh = aderor/gamold
        adgamold = adgamold-aderor*(abs(gamold-gamaan)/(gamold*gamold))
        adgamaan = adgamaan-aderorh*sign(1.d0,gamold-gamaan)
        adgamold = adgamold+aderorh*sign(1.d0,gamold-gamaan)
        aderor = 0.
        adgams(2,2) = adgams(2,2)+adgamaan
        adgamaan = 0.
        call adactcof( cat,an,adcat,adan,adgams )
        adan(3) = 0.
        adman = adman+adan(2)
        adan(2) = 0.
        admas = admas+adan(1)
        adan(1) = 0.
        admnh4 = admnh4+adcat(2)
        adcat(2) = 0.
        adcat(1) = 0.
        admnh4 = admnh4+adynh4*wh2o
        adwh2o = adwh2o+adynh4*mnh4
        adynh4 = 0.
        adman = adman+admnh4
        admas = admas+2*admnh4
        admnh4 = 0.
        adtso4 = adtso4+admas/wh2o
        adwh2o = adwh2o-admas*(tso4/(wh2o*wh2o))
        admas = 0.
        adwh2o = adwh2o-adman*(xno3/(wh2o*wh2o))
        adxno3 = adxno3+adman/wh2o
        adman = 0.
        adah2o = adah2o+0.001*adwh2o
        adwh2o = 0.
        ynh4 = ynh4h
        call adawater( irh,tso4,ynh4,xno3,adtso4,adynh4,adxno3,adah2o )
        if (aa .ne. 0.d0) then
          dd = sqrt(disc)
          xxq = -(0.5d0*(bb+sign(1.d0,bb)*dd))
          rr1 = xxq/aa
          rr2 = cc/xxq
          if (rr1*rr2 .lt. 0.d0) then
            xno3 = max(rr1,rr2)
          else
            xno3 = min(rr1,rr2)
          endif
        else
          xno3 = -(cc/bb)
        endif
        adtno3 = adtno3+adxno3*(0.5-sign(0.5d0,tno3-xno3))
        adxno3 = adxno3*(0.5+sign(0.5d0,tno3-xno3))
        if (aa .ne. 0.d0) then
          if (rr1*rr2 .lt. 0.d0) then
            adrr1 = adrr1+adxno3*(0.5+sign(0.5d0,rr1-rr2))
            adrr2 = adrr2+adxno3*(0.5-sign(0.5d0,rr1-rr2))
            adxno3 = 0.
          else
            adrr1 = adrr1+adxno3*(0.5+sign(0.5d0,rr2-rr1))
            adrr2 = adrr2+adxno3*(0.5-sign(0.5d0,rr2-rr1))
            adxno3 = 0.
          endif
          adcc = adcc+adrr2/xxq
          adxxq = adxxq-adrr2*(cc/(xxq*xxq))
          adrr2 = 0.
          adaa = adaa-adrr1*(xxq/(aa*aa))
          adxxq = adxxq+adrr1/aa
          adrr1 = 0.
          adbb = adbb-0.5d0*adxxq
          addd = addd-0.5d0*adxxq*sign(1.d0,bb)
          adxxq = 0.
          addisc = addisc+addd*(1./(2.*sqrt(disc)))
          addd = 0.
        else
          adbb = adbb+adxno3*(cc/(bb*bb))
          adcc = adcc-adxno3/bb
          adxno3 = 0.
        endif
        adaa = adaa-4*addisc*cc
        adbb = adbb+2*addisc*bb
        adcc = adcc-4*addisc*aa
        addisc = 0.
        adkw2 = adkw2-adcc*tno3*(tnh4-twoso4)
        adtnh4 = adtnh4-adcc*kw2*tno3
        adtno3 = adtno3-adcc*kw2*(tnh4-twoso4)
        adtwoso4 = adtwoso4+adcc*kw2*tno3
        adcc = 0.
        adkw2 = adkw2+adbb*(tno3+tnh4-twoso4)
        adtnh4 = adtnh4+adbb*kw2
        adtno3 = adtno3+adbb*kw2
        adtwoso4 = adtwoso4+adbb*(1-kw2)
        adbb = 0.
        adkw2 = adkw2-adaa
        adaa = 0.
        adgasqd = adgasqd-adkw2*(kan*wsqd/(gasqd*gasqd))
        adwsqd = adwsqd+adkw2*(kan/gasqd)
        adkw2 = 0.
        wh2o = wh2oh
        adwh2o = adwh2o+2*adwsqd*wh2o
        adwsqd = 0.
        gamaan = gamaanh
        adgamaan = adgamaan+2*adgasqd*gamaan
        adgasqd = 0.
      end do
      adtwoso4 = adtwoso4+adynh4
      adynh4 = 0.
      adynh4 = adynh4+adanh4*mwnh4
      adanh4 = 0.
      adano3 = 0.
      adtso4 = adtso4+adaso4*mwso4
      adaso4 = 0.
      adah2o = adah2o+0.001d0*adwh2o
      adwh2o = 0.
      ynh4 = twoso4
      call adawater( irh,tso4,ynh4,tno3,adtso4,adynh4,adtno3,adah2o )
      adtwoso4 = adtwoso4+adynh4
      adynh4 = 0.
      adtso4 = adtso4+2*adtwoso4
      adtwoso4 = 0.
      ano3 = in(4)
      adano3 = adano3+adtmasshno3*(0.5-sign(0.5d0,0.d0-(gno3+ano3)))
      adgno3 = adgno3+adtmasshno3*(0.5-sign(0.5d0,0.d0-(gno3+ano3)))
      adtmasshno3 = 0.
      adanh4 = adanh4+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh4)
      adgnh3 = adgnh3+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh3)
      adtnh4 = 0.
      adano3 = adano3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwno3)
      adgno3 = adgno3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwhno3)
      adtno3 = 0.
      adso4 = adso4+adtso4*((0.5-sign(0.5d0,floor-so4/mwso4))/mwso4)
      adtso4 = 0.
      adin(5) = adin(5)+adanh4
      adanh4 = 0.
      adin(4) = adin(4)+adano3
      adano3 = 0.
      adin(3) = adin(3)+adgnh3
      adgnh3 = 0.
      adin(2) = adin(2)+adgno3
      adgno3 = 0.
      adin(1) = adin(1)+adso4
      adso4 = 0.

      end SUBROUTINE ADRPMARES_6

!-----------------------------------------------------------------------------


C                           DISCLAIMER
C
C   This file was generated by TAMC version 5.3.2
C
C   THE AUTHOR DOES NOT MAKE  ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR  RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
C   OR USEFULNESS  OF ANY INFORMATION OR PROCESS  DISCLOSED, OR REPRESENTS
C   THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
C   THIS CODE IS FOR NON-PROFIT-ORIENTATED ACADEMIC RESEARCH AND EDUCATION
C   ONLY.  ANY COMMERCIAL OR OTHER PROFIT-ORIENTATED USE OR  EVALUATION IS
C   STRICTLY  FORBIDDEN.  PASSING  THIS CODE  TO  ANY THIRD  PARTY IS  NOT
C   ALLOWED.
C
C   FOR COMMERCIAL OR  OTHER PROFIT-ORIENTATED APPLICATIONS PLEASE CONTACT
C   info@FastOpt.com
C
      subroutine adrpmares_7( in, par, adin, adout, 
     &                         I,  J,   L )

C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================

C==============================================
C define parameters                            
C==============================================
      TYPE (XPLEX) floor
      parameter ( floor = xplex(1.d-30,0d0) )
      TYPE (XPLEX) mwhno3
      parameter ( mwhno3 = xplex(63.01287d0,0d0) )
      TYPE (XPLEX) mwnh3
      parameter ( mwnh3 = xplex(17.03061d0,0d0) )
      TYPE (XPLEX) mwnh4
      parameter ( mwnh4 = xplex(18.03858d0,0d0) )
      TYPE (XPLEX) mwno3
      parameter ( mwno3 = xplex(62.0049d0,0d0) )
      TYPE (XPLEX) mwso4
      parameter ( mwso4 = xplex(96.0576d0,0d0) )

C==============================================
C define common blocks
C==============================================
C==============================================
C define arguments
C==============================================
      TYPE (XPLEX) adin(5)
      TYPE (XPLEX) adout(8)
      TYPE (XPLEX) in(5)
      TYPE (XPLEX) par(2)
      INTEGER :: I, J, L

C==============================================
C define local variables
C==============================================
      TYPE (XPLEX) adah2o
      TYPE (XPLEX) adanh4
      TYPE (XPLEX) adano3
      TYPE (XPLEX) adano3_in
      TYPE (XPLEX) adaso4
      TYPE (XPLEX) adgnh3
      TYPE (XPLEX) adgno3
      TYPE (XPLEX) adgno3_in
      TYPE (XPLEX) adso4
      TYPE (XPLEX) adtnh4
      TYPE (XPLEX) adtno3
      TYPE (XPLEX) adtso4
      TYPE (XPLEX) adtwoso4
      TYPE (XPLEX) adxno3
      TYPE (XPLEX) adynh4
      TYPE (XPLEX) anh4
      TYPE (XPLEX) ano3
      TYPE (XPLEX) gnh3
      TYPE (XPLEX) gno3
      integer irh
      TYPE (XPLEX) rh
      TYPE (XPLEX) so4
      TYPE (XPLEX) tnh4
      TYPE (XPLEX) tno3
      TYPE (XPLEX) tso4
      TYPE (XPLEX) twoso4
      TYPE (XPLEX) xno3
      TYPE (XPLEX) ynh4

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adah2o = 0.
      adanh4 = 0.
      adano3 = 0.
      adano3_in = 0.
      adaso4 = 0.
      adgnh3 = 0.
      adgno3 = 0.
      adgno3_in = 0.
      adso4 = 0.
      adtnh4 = 0.
      adtno3 = 0.
      adtso4 = 0.
      adtwoso4 = 0.
      adxno3 = 0.
      adynh4 = 0.

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      so4 = in(1)
      gno3 = in(2)
      gnh3 = in(3)
      ano3 = in(4)
      anh4 = in(5)
      rh = par(1)
      tso4 = max(floor,so4/mwso4)
      tno3 = max(0.d0,ano3/mwno3+gno3/mwhno3)
      tnh4 = max(0.d0,gnh3/mwnh3+anh4/mwnh4)
      irh = nint(100.*rh)
      irh = max(1,irh)
      irh = min(99,irh)
      twoso4 = 2.*tso4
      xno3 = tno3/mwno3
      ynh4 = twoso4
      adgnh3 = adgnh3+adout(8)
      adout(8) = 0.d0
      adgno3 = adgno3+adout(7)
      adout(7) = 0.d0
      adso4 = adso4+adout(6)
      adout(6) = 0.d0
      adanh4 = adanh4+adout(5)
      adout(5) = 0.d0
      adah2o = adah2o+adout(4)
      adout(4) = 0.d0
      adano3 = adano3+adout(3)
      adout(3) = 0.d0
      adout(2) = 0.d0
      adaso4 = adaso4+adout(1)
      adout(1) = 0.d0
      adtnh4 = adtnh4+adgnh3*(0.5-sign(0.5d0,floor-mwnh3*(tnh4-ynh4)))*
     $mwnh3
      adynh4 = adynh4-adgnh3*(0.5-sign(0.5d0,floor-mwnh3*(tnh4-ynh4)))*
     $mwnh3
      adgnh3 = 0.
      adano3_in = adano3_in+adano3
      adano3 = 0.
      adgno3_in = adgno3_in+adgno3
      adgno3 = 0.
      call adawater( irh,tso4,ynh4,xno3,adtso4,adynh4,adxno3,adah2o )
      adynh4 = adynh4+adanh4*mwnh4
      adanh4 = 0.
      adtwoso4 = adtwoso4+adynh4
      adynh4 = 0.
      adtno3 = adtno3+adxno3/mwno3
      adxno3 = 0.
      adtso4 = adtso4+adaso4*mwso4
      adaso4 = 0.
      ynh4 = twoso4
      call adawater( irh,tso4,ynh4,tno3,adtso4,adynh4,adtno3,adah2o )
      adtwoso4 = adtwoso4+adynh4
      adynh4 = 0.
      adtso4 = adtso4+2*adtwoso4
      adtwoso4 = 0.
      adano3 = adano3+adano3_in
      adano3_in = 0.
      adgno3 = adgno3+adgno3_in
      adgno3_in = 0.
      adanh4 = adanh4+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh4)
      adgnh3 = adgnh3+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh3)
      adtnh4 = 0.
      adano3 = adano3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwno3)
      adgno3 = adgno3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwhno3)
      adtno3 = 0.
      adso4 = adso4+adtso4*((0.5-sign(0.5d0,floor-so4/mwso4))/mwso4)
      adtso4 = 0.
      adin(5) = adin(5)+adanh4
      adanh4 = 0.
      adin(4) = adin(4)+adano3
      adano3 = 0.
      adin(3) = adin(3)+adgnh3
      adgnh3 = 0.
      adin(2) = adin(2)+adgno3
      adgno3 = 0.
      adin(1) = adin(1)+adso4
      adso4 = 0.

      end SUBROUTINE ADRPMARES_7

!-----------------------------------------------------------------------------


C                           DISCLAIMER
C
C   This file was generated by TAMC version 5.3.2
C
C   THE AUTHOR DOES NOT MAKE  ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR  RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
C   OR USEFULNESS  OF ANY INFORMATION OR PROCESS  DISCLOSED, OR REPRESENTS
C   THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
C   THIS CODE IS FOR NON-PROFIT-ORIENTATED ACADEMIC RESEARCH AND EDUCATION
C   ONLY.  ANY COMMERCIAL OR OTHER PROFIT-ORIENTATED USE OR  EVALUATION IS
C   STRICTLY  FORBIDDEN.  PASSING  THIS CODE  TO  ANY THIRD  PARTY IS  NOT
C   ALLOWED.
C
C   FOR COMMERCIAL OR  OTHER PROFIT-ORIENTATED APPLICATIONS PLEASE CONTACT
C   info@FastOpt.com
C
      subroutine adrpmares_8( in, par, adin, adout,
     &                         I,  J,   L )

C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================

C==============================================
C define parameters                            
C==============================================
      TYPE (XPLEX) floor
      parameter ( floor = xplex(1.d-30,0d0) )
      TYPE (XPLEX) mwhno3
      parameter ( mwhno3 = xplex(63.01287d0,0d0) )
      TYPE (XPLEX) mwnh3
      parameter ( mwnh3 = xplex(17.03061d0,0d0) )
      TYPE (XPLEX) mwnh4
      parameter ( mwnh4 = xplex(18.03858d0,0d0) )
      TYPE (XPLEX) mwno3
      parameter ( mwno3 = xplex(62.0049d0,0d0) )
      TYPE (XPLEX) mwso4
      parameter ( mwso4 = xplex(96.0576d0,0d0) )

C==============================================
C define common blocks
C==============================================
C==============================================
C define arguments
C==============================================
      TYPE (XPLEX) adin(5)
      TYPE (XPLEX) adout(8)
      TYPE (XPLEX) in(5)
      TYPE (XPLEX) par(2)
      INTEGER :: I, J, L

C==============================================
C define local variables
C==============================================
      TYPE (XPLEX) adah2o
      TYPE (XPLEX) adahso4
      TYPE (XPLEX) adanh4
      TYPE (XPLEX) adano3
      TYPE (XPLEX) adano3_in
      TYPE (XPLEX) adaso4
      TYPE (XPLEX) adgnh3
      TYPE (XPLEX) adgno3
      TYPE (XPLEX) adso4
      TYPE (XPLEX) adtmasshno3
      TYPE (XPLEX) adtnh4
      TYPE (XPLEX) adtno3
      TYPE (XPLEX) adtso4
      TYPE (XPLEX) anh4
      TYPE (XPLEX) ano3
      TYPE (XPLEX) gnh3
      TYPE (XPLEX) gno3
      integer irh
      TYPE (XPLEX) rh
      TYPE (XPLEX) so4
      TYPE (XPLEX) tnh4
      TYPE (XPLEX) tno3
      TYPE (XPLEX) tso4

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adah2o = 0.
      adahso4 = 0.
      adanh4 = 0.
      adano3 = 0.
      adano3_in = 0.
      adaso4 = 0.
      adgnh3 = 0.
      adgno3 = 0.
      adso4 = 0.
      adtmasshno3 = 0.
      adtnh4 = 0.
      adtno3 = 0.
      adtso4 = 0.

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      so4 = in(1)
      gno3 = in(2)
      gnh3 = in(3)
      ano3 = in(4)
      anh4 = in(5)
      rh = par(1)
      tso4 = max(floor,so4/mwso4)
      tno3 = max(0.d0,ano3/mwno3+gno3/mwhno3)
      tnh4 = max(0.d0,gnh3/mwnh3+anh4/mwnh4)
      irh = nint(100.*rh)
      irh = max(1,irh)
      irh = min(99,irh)
      adgnh3 = adgnh3+adout(8)
      adout(8) = 0.d0
      adgno3 = adgno3+adout(7)
      adout(7) = 0.d0
      adso4 = adso4+adout(6)
      adout(6) = 0.d0
      adanh4 = adanh4+adout(5)
      adout(5) = 0.d0
      adah2o = adah2o+adout(4)
      adout(4) = 0.d0
      adano3 = adano3+adout(3)
      adout(3) = 0.d0
      adahso4 = adahso4+adout(2)
      adout(2) = 0.d0
      adaso4 = adaso4+adout(1)
      adout(1) = 0.d0
      adgnh3 = 0.
      adano3 = adano3-adgno3
      adtmasshno3 = adtmasshno3+adgno3
      adgno3 = 0.
      adano3_in = adano3_in+adano3
      adano3 = 0.
      adtnh4 = adtnh4+adanh4*mwnh4
      adanh4 = 0.
      adtso4 = adtso4+adahso4*mwso4
      adahso4 = 0.
      adaso4 = 0.
      call adawater( irh,tso4,tnh4,tno3,adtso4,adtnh4,adtno3,adah2o )
      adano3 = adano3+adtmasshno3*(0.5-sign(0.5d0,0.d0-(gno3+ano3)))
      adgno3 = adgno3+adtmasshno3*(0.5-sign(0.5d0,0.d0-(gno3+ano3)))
      adtmasshno3 = 0.
      adano3 = adano3+adano3_in
      adano3_in = 0.
      adanh4 = adanh4+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh4)
      adgnh3 = adgnh3+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh3)
      adtnh4 = 0.
      adano3 = adano3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwno3)
      adgno3 = adgno3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwhno3)
      adtno3 = 0.
      adso4 = adso4+adtso4*((0.5-sign(0.5d0,floor-so4/mwso4))/mwso4)
      adtso4 = 0.
      adin(5) = adin(5)+adanh4
      adanh4 = 0.
      adin(4) = adin(4)+adano3
      adano3 = 0.
      adin(3) = adin(3)+adgnh3
      adgnh3 = 0.
      adin(2) = adin(2)+adgno3
      adgno3 = 0.
      adin(1) = adin(1)+adso4
      adso4 = 0.

      end SUBROUTINE ADRPMARES_8
!------------------------------------------------------------------------------

      subroutine adcubic( a2, a1, a0, ada2, ada1, ada0, adcrutes )
C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================
      USE MYTYPE
      USE COMPLEXIFY
      implicit none

C==============================================
C define parameters
C==============================================
      TYPE (XPLEX) one
      parameter ( one = xplex(1.d0,0d0) )
      TYPE (XPLEX) one3rd
      parameter ( one3rd = xplex(0.333333333d0,0d0) )
      TYPE (XPLEX) sqrt3
      parameter ( sqrt3 = xplex(1.732050808d0,0d0) )

C==============================================
C define common blocks
C==============================================
C==============================================
C define arguments
C==============================================
      TYPE (XPLEX) a0
      TYPE (XPLEX) a1
      TYPE (XPLEX) a2
      TYPE (XPLEX) ada0
      TYPE (XPLEX) ada1
      TYPE (XPLEX) ada2
      TYPE (XPLEX) adcrutes(3)

C==============================================
C define local variables
C==============================================
      TYPE (XPLEX) a2sq
      TYPE (XPLEX) ada2sq
      TYPE (XPLEX) adcosth
      TYPE (XPLEX) addum1
      TYPE (XPLEX) addum2
      TYPE (XPLEX) adpart1
      TYPE (XPLEX) adpart2
      TYPE (XPLEX) adpart3
      TYPE (XPLEX) adphi
      TYPE (XPLEX) adqq
      TYPE (XPLEX) adrr
      TYPE (XPLEX) adrrsq
      TYPE (XPLEX) adsinth
      TYPE (XPLEX) adtheta
      TYPE (XPLEX) adyy1
      TYPE (XPLEX) adyy2
      TYPE (XPLEX) adyy3
      TYPE (XPLEX) costh
      TYPE (XPLEX) crutes(3)
      TYPE (XPLEX) dum1
      TYPE (XPLEX) dum2
      TYPE (XPLEX) part1
      TYPE (XPLEX) part2
      TYPE (XPLEX) part3
      TYPE (XPLEX) phi
      TYPE (XPLEX) qq
      TYPE (XPLEX) rr
      TYPE (XPLEX) rrsq
      TYPE (XPLEX) sinth
      TYPE (XPLEX) theta
      TYPE (XPLEX) yy1
      TYPE (XPLEX) yy2
      TYPE (XPLEX) yy3

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      ada2sq = 0.
      adcosth = 0.
      addum1 = 0.
      addum2 = 0.
      adpart1 = 0.
      adpart2 = 0.
      adpart3 = 0.
      adphi = 0.
      adqq = 0.
      adrr = 0.
      adrrsq = 0.
      adsinth = 0.
      adtheta = 0.
      adyy1 = 0.
      adyy2 = 0.
      adyy3 = 0.

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      a2sq = a2*a2
      qq = (a2sq-3.d0*a1)/9.d0
      rr = (a2*(2.d0*a2sq-9.d0*a1)+27.d0*a0)/54.d0
      dum1 = qq*qq*qq
      rrsq = rr*rr
      dum2 = dum1-rrsq
      if (dum2 .ge. 0.d0) then
        phi = sqrt(dum1)
        if (abs(phi) .lt. 1.d-20) then
          crutes(1) = 0.d0
          crutes(2) = 0.d0
          crutes(3) = 0.d0
        endif
        theta = acos(rr/phi)/3.d0
        costh = cos(theta)
        sinth = sin(theta)
        part1 = sqrt(qq)
        yy1 = part1*costh
        yy2 = yy1-a2/3.d0
        yy3 = sqrt3*part1*sinth
        crutes(3) = (-(2.d0*yy1))-a2/3.d0
        crutes(2) = yy2+yy3
        crutes(1) = yy2-yy3
        if (crutes(1) .lt. 0.d0) then
          crutes(1) = 1.d+9
        endif
        if (crutes(2) .lt. 0.d0) then
          crutes(2) = 1.d+9
        endif
        if (crutes(3) .lt. 0.d0) then
          crutes(3) = 1.d+9
        endif
        adcrutes(2) = adcrutes(2)+adcrutes(1)*(0.5-sign(0.5d0,crutes(2)-
     $crutes(1)))
        adcrutes(1) = adcrutes(1)*(0.5+sign(0.5d0,crutes(2)-crutes(1)))
        crutes(1) = yy2-yy3
        if (crutes(1) .lt. 0.d0) then
          crutes(1) = 1.d+9
        endif
        if (crutes(2) .lt. 0.d0) then
          crutes(2) = 1.d+9
        endif
        if (crutes(3) .lt. 0.d0) then
          adcrutes(3) = 0.
        endif
        crutes(1) = yy2-yy3
        if (crutes(1) .lt. 0.d0) then
          crutes(1) = 1.d+9
        endif
        if (crutes(2) .lt. 0.d0) then
          adcrutes(2) = 0.
        endif
        crutes(1) = yy2-yy3
        if (crutes(1) .lt. 0.d0) then
          adcrutes(1) = 0.
        endif
        adyy2 = adyy2+adcrutes(1)
        adyy3 = adyy3-adcrutes(1)
        adcrutes(1) = 0.
        adyy2 = adyy2+adcrutes(2)
        adyy3 = adyy3+adcrutes(2)
        adcrutes(2) = 0.
        ada2 = ada2-0.333333333333d0*adcrutes(3)
        adyy1 = adyy1-2*adcrutes(3)
        adcrutes(3) = 0.
        adpart1 = adpart1+adyy3*sqrt3*sinth
        adsinth = adsinth+adyy3*sqrt3*part1
        adyy3 = 0.
        ada2 = ada2-0.333333333333d0*adyy2
        adyy1 = adyy1+adyy2
        adyy2 = 0.
        adcosth = adcosth+adyy1*part1
        adpart1 = adpart1+adyy1*costh
        adyy1 = 0.
        adqq = adqq+adpart1*(1./(2.*sqrt(qq)))
        adpart1 = 0.
        adtheta = adtheta+adsinth*cos(theta)
        adsinth = 0.
        adtheta = adtheta-adcosth*sin(theta)
        adcosth = 0.
        adphi = adphi+adtheta*(1./sqrt(1.-(rr/phi)**2)*(rr/(phi*phi))/
     $3.d0)
        adrr = adrr-adtheta*(1./sqrt(1.-(rr/phi)**2)/phi/3.d0)
        adtheta = 0.
        if (abs(phi) .lt. 1.d-20) then
          adcrutes(3) = 0.
          adcrutes(2) = 0.
          adcrutes(1) = 0.
        endif
        addum1 = addum1+adphi*(1./(2.*sqrt(dum1)))
        adphi = 0.
      else
        part1 = sqrt(rrsq-dum1)
        part2 = abs(rr)
        part3 = (part1+part2)**one3rd
        adcrutes(3) = 0.
        adcrutes(2) = 0.
        ada2 = ada2-0.333333333333d0*adcrutes(1)
        adpart3 = adpart3-adcrutes(1)*(1-qq/(part3*part3))*sign(one,rr)
        adqq = adqq-adcrutes(1)/part3*sign(one,rr)
        adcrutes(1) = 0.
        adpart1 = adpart1+adpart3*one3rd*(part1+part2)**(one3rd-1)
        adpart2 = adpart2+adpart3*one3rd*(part1+part2)**(one3rd-1)
        adpart3 = 0.
        adrr = adrr+adpart2*sign(1.d0,rr)
        adpart2 = 0.
        addum1 = addum1-adpart1*(1./(2.*sqrt(rrsq-dum1)))
        adrrsq = adrrsq+adpart1*(1./(2.*sqrt(rrsq-dum1)))
        adpart1 = 0.
      endif
      addum1 = addum1+addum2
      adrrsq = adrrsq-addum2
      addum2 = 0.
      adrr = adrr+2*adrrsq*rr
      adrrsq = 0.
      adqq = adqq+3*addum1*qq*qq
      addum1 = 0.
      ada0 = ada0+0.5d0*adrr
      ada1 = ada1+adrr*((-9)*a2/54.d0)
      ada2 = ada2+adrr*((2.d0*a2sq-9.d0*a1)/54.d0)
      ada2sq = ada2sq+adrr*(2*a2/54.d0)
      adrr = 0.
      ada1 = ada1-0.333333333333d0*adqq
      ada2sq = ada2sq+0.111111111111d0*adqq
      adqq = 0.
      ada2 = ada2+2*ada2sq*a2
      ada2sq = 0.

      end SUBROUTINE ADCUBIC 
!------------------------------------------------------------------------------

      subroutine adrpmares_6_D5( in, par, adin, adout,
     &                      I, J, L )

!
!******************************************************************************
!  Subroutine adrpmares_6_D5 was created using a modified version of 
!  rpmares_short6.f which replaced the RETURN structure with a DOWHILE loop. 
!  (dkh, 06/01/05)
! 
! Notes
! (1 ) The following changes are made to the code returned by TAMC:
!       - Change the routine name
!       - Expand argument list to include I, J, L
!       - Eliminate the OUT variable as we are not returning results of fwd 
!          calculation from this subroutine. 
!       - Replace reference to TAMC storage routines with reference to our 
!          checkpointing variables xxx_fwd which are initialized in 
!          RECOMP_RPMARES.  Comment out the portions of the code that were used
!          to recompute these variables (idow, gamaan, wh2o, ynh4).
!  
! (2 ) Unlike previous version, don't bother to use checkpointed values of 
!       gamanold Can just use gamaan_fwd(I,J,L,idow-1) to restore gamanold.
!
!******************************************************************************
!

      ! Reference to f90 modules
      USE CHECKPT_MOD
      USE RPMARES_MOD,   ONLY : CUBIC, AWATER, ACTCOF

#     include "CMN_SIZE"  ! Size params

C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================
!      USE MYTYPE
!      USE COMPLEXIFY
!      implicit none

C==============================================
C define parameters                            
C==============================================
      TYPE (XPLEX) floor
      parameter ( floor = xplex(1.d-30,0d0) )
      TYPE (XPLEX) mwhno3
      parameter ( mwhno3 = xplex(63.01287d0,0d0) )
      TYPE (XPLEX) mwnh3
      parameter ( mwnh3 = xplex(17.03061d0,0d0) )
      TYPE (XPLEX) mwnh4
      parameter ( mwnh4 = xplex(18.03858d0,0d0) )
      TYPE (XPLEX) mwno3
      parameter ( mwno3 = xplex(62.0049d0,0d0) )
      TYPE (XPLEX) mwso4
      parameter ( mwso4 = xplex(96.0576d0,0d0) )

C==============================================
C define common blocks
C==============================================
C==============================================
C define arguments
C==============================================
      TYPE (XPLEX) adin(5)
      TYPE (XPLEX) adout(8)
      TYPE (XPLEX) in(5)
!      TYPE (XPLEX) out(8)
      TYPE (XPLEX) par(2)

C==============================================
C define local variables
C==============================================
      TYPE (XPLEX) aa
      TYPE (XPLEX) adaa
      TYPE (XPLEX) adah2o
      TYPE (XPLEX) adan(3)
      TYPE (XPLEX) adanh4
      TYPE (XPLEX) adano3
      TYPE (XPLEX) adaso4
      TYPE (XPLEX) adbb
      TYPE (XPLEX) adcat(2)
      TYPE (XPLEX) adcc
      TYPE (XPLEX) addd
      TYPE (XPLEX) addisc
      TYPE (XPLEX) aderor
      TYPE (XPLEX) aderori
      TYPE (XPLEX) adgamaan
      TYPE (XPLEX) adgamold
      TYPE (XPLEX) adgams(2,3)
      TYPE (XPLEX) adgasqd
      TYPE (XPLEX) adgnh3
      TYPE (XPLEX) adgno3
      TYPE (XPLEX) adkw2
      TYPE (XPLEX) adman
      TYPE (XPLEX) admas
      TYPE (XPLEX) admnh4
      TYPE (XPLEX) adrr1
      TYPE (XPLEX) adrr2
      TYPE (XPLEX) adso4
      TYPE (XPLEX) adtmasshno3
      TYPE (XPLEX) adtnh4
      TYPE (XPLEX) adtno3
      TYPE (XPLEX) adtso4
      TYPE (XPLEX) adtwoso4
      TYPE (XPLEX) adwh2o
      TYPE (XPLEX) adwsqd
      TYPE (XPLEX) adxno3
      TYPE (XPLEX) adxxq
      TYPE (XPLEX) adynh4
      TYPE (XPLEX) ah2o
      TYPE (XPLEX) ahso4
      TYPE (XPLEX) an(3)
      TYPE (XPLEX) anh4
      TYPE (XPLEX) ano3
      TYPE (XPLEX) aso4
      TYPE (XPLEX) bb
      TYPE (XPLEX) cat(2)
      TYPE (XPLEX) cc
      logical converged
      TYPE (XPLEX) dd
      TYPE (XPLEX) disc
      TYPE (XPLEX) eror
      TYPE (XPLEX) gamaan
      TYPE (XPLEX) gamaani
      TYPE (XPLEX) gamold
      TYPE (XPLEX) gams(2,3)
      TYPE (XPLEX) gasqd
      TYPE (XPLEX) gnh3
      TYPE (XPLEX) gno3
      integer idow
      integer idow2
      integer ip1
      integer ip2
      integer irh
      TYPE (XPLEX) k1a
      TYPE (XPLEX) kan
      TYPE (XPLEX) khat
      TYPE (XPLEX) kna
      TYPE (XPLEX) kph
      TYPE (XPLEX) kw
      TYPE (XPLEX) kw2
      TYPE (XPLEX) man
      TYPE (XPLEX) mas
      TYPE (XPLEX) mnh4
      TYPE (XPLEX) molnu
      integer ndow
      integer nnn
      TYPE (XPLEX) phibar
      TYPE (XPLEX) rh
      TYPE (XPLEX) rr1
      TYPE (XPLEX) rr2
      TYPE (XPLEX) so4
      TYPE (XPLEX) t1
      TYPE (XPLEX) t2
      TYPE (XPLEX) t3
      TYPE (XPLEX) t4
      TYPE (XPLEX) t6
      TYPE (XPLEX) temp
      TYPE (XPLEX) tmasshno3
      TYPE (XPLEX) tnh4
      TYPE (XPLEX) tno3
      TYPE (XPLEX) toler1
      TYPE (XPLEX) tso4
      TYPE (XPLEX) twoso4
      TYPE (XPLEX) wh2o
      TYPE (XPLEX) wh2oi
      TYPE (XPLEX) wsqd
      TYPE (XPLEX) xno3
      TYPE (XPLEX) xxq
      TYPE (XPLEX) ynh4
      TYPE (XPLEX) ynh4i
      INTEGER I, J, L

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adaa = 0.
      adah2o = 0.
      do ip1 = 1, 3
        adan(ip1) = 0.
      end do
      adanh4 = 0.
      adano3 = 0.
      adaso4 = 0.
      adbb = 0.
      do ip1 = 1, 2
        adcat(ip1) = 0.
      end do
      adcc = 0.
      addd = 0.
      addisc = 0.
      aderor = 0.
      adgamaan = 0.
      adgamold = 0.
      do ip2 = 1, 3
        do ip1 = 1, 2
          adgams(ip1,ip2) = 0.
        end do
      end do
      adgasqd = 0.
      adgnh3 = 0.
      adgno3 = 0.
      adkw2 = 0.
      adman = 0.
      admas = 0.
      admnh4 = 0.
      adrr1 = 0.
      adrr2 = 0.
      adso4 = 0.
      adtmasshno3 = 0.
      adtnh4 = 0.
      adtno3 = 0.
      adtso4 = 0.
      adtwoso4 = 0.
      adwh2o = 0.
      adwsqd = 0.
      adxno3 = 0.
      adxxq = 0.
      adynh4 = 0.

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
C----------------------------------------------
C FUNCTION AND TAPE COMPUTATIONS
C----------------------------------------------
      so4 = in(1)
      gno3 = in(2)
      gnh3 = in(3)
      ano3 = in(4)
      anh4 = in(5)
      rh = par(1)
      temp = par(2)
      tso4 = max(floor,so4/mwso4)
      tno3 = max(0.d0,ano3/mwno3+gno3/mwhno3)
      tnh4 = max(0.d0,gnh3/mwnh3+anh4/mwnh4)
      tmasshno3 = max(0.d0,gno3+ano3)
      irh = nint(100.*rh)
      irh = max(1,irh)
      irh = min(99,irh)
      t6 = 8.2d-11*temp
      t1 = 298.d0/temp
      t2 = log(t1)
      t3 = t1-1.d0
      t4 = 1.d0+t2-t1
      kna = 2511000.d0*exp(29.17d0*t3+16.83d0*t4)*t6
      k1a = 0.00001805d0*exp((-(1.5d0*t3))+26.92d0*t4)
      kw = 1.01d-14*exp((-(22.52d0*t3))+26.92d0*t4)
      kph = 57.639d0*exp(13.79d0*t3-5.39d0*t4)*t6
      khat = kph*k1a/kw
      kan = kna*khat
      toler1 = 0.00001d0
      gamold = 1.d0
      gamaan = 0.1
      twoso4 = 2.*tso4
      ynh4 = twoso4
      call awater( irh,tso4,ynh4,tno3,ah2o )
      wh2o = 0.001d0*ah2o
      aso4 = tso4*mwso4
      ahso4 = 0.d0
      ano3 = 0.d0
      anh4 = ynh4*mwnh4
      ynh4 = twoso4
      nnn = 0
      converged =  .false. 
      idow = 0
      !The following section is used by TAMC to recompute idow. 
      !We comment this out and use the values from  nitr_max
!      do while (converged .eq.  .false.  .or. nnn .gt. 50 )
!        idow = idow+1
!        nnn = nnn+1
!        gasqd = gamaan*gamaan
!        wsqd = wh2o*wh2o
!        kw2 = kan*wsqd/gasqd
!        aa = 1.-kw2
!        bb = twoso4+kw2*(tno3+tnh4-twoso4)
!        cc = -(kw2*tno3*(tnh4-twoso4))
!        disc = bb*bb-4.*aa*cc
!        if (aa .ne. 0.d0) then
!          dd = sqrt(disc)
!          xxq = -(0.5d0*(bb+sign(1.d0,bb)*dd))
!          rr1 = xxq/aa
!          rr2 = cc/xxq
!          if (rr1*rr2 .lt. 0.d0) then
!            xno3 = max(rr1,rr2)
!          else
!            xno3 = min(rr1,rr2)
!          endif
!        else
!          xno3 = -(cc/bb)
!        endif
!        xno3 = min(xno3,tno3)
!        call awater( irh,tso4,ynh4,xno3,ah2o )
!        wh2o = 0.001*ah2o
!        man = xno3/wh2o
!        mas = tso4/wh2o
!        mnh4 = 2.*mas+man
!        ynh4 = mnh4*wh2o
!        cat(1) = 0.
!        cat(2) = mnh4
!        an(1) = mas
!        an(2) = man
!        an(3) = 0.
!        call actcof( cat,an,gams,molnu,phibar )
!        gamaan = gams(2,2)
!        eror = abs(gamold-gamaan)/gamold
!        gamold = gamaan
!        if (eror .le. toler1) then
!          aso4 = tso4*mwso4
!          ahso4 = 0.
!          ano3 = xno3*mwno3
!          anh4 = ynh4*mwnh4
!          gno3 = max(floor,tmasshno3-ano3)
!          gnh3 = mwnh3*max(floor,tnh4-ynh4)
!          converged =  .true. 
!        endif
!      end do
!      call adstore( 'memory_1_rpmares_idow',21,idow,4,1,1 )
!      out(1) = aso4
!      out(2) = ahso4
!      out(3) = ano3
!      out(4) = ah2o
!      out(5) = anh4
!      out(6) = so4
!      out(7) = gno3
!      out(8) = gnh3

      idow = nitr_max(I,J,L)

C----------------------------------------------
C ADJOINT COMPUTATIONS
C----------------------------------------------
      ano3 = in(4)
      adgnh3 = adgnh3+adout(8)
      adout(8) = 0.d0
      adgno3 = adgno3+adout(7)
      adout(7) = 0.d0
      adso4 = adso4+adout(6)
      adout(6) = 0.d0
      adanh4 = adanh4+adout(5)
      adout(5) = 0.d0
      adah2o = adah2o+adout(4)
      adout(4) = 0.d0
      adano3 = adano3+adout(3)
      adout(3) = 0.d0
      adout(2) = 0.d0
      adaso4 = adaso4+adout(1)
      adout(1) = 0.d0
      !call adresto( 'memory_1_rpmares_idow',21,idow,4,1,1 )
      ndow = idow
      do idow = ndow, 1, -1
        gamaan = 0.1
        ynh4 = twoso4
        call awater( irh,tso4,ynh4,tno3,ah2o )
        wh2o = 0.001d0*ah2o
        ynh4 = twoso4
         ! The following section is used to recompute gamaan, ynh4 and 
         ! wh20 at each iteration.  Instead, comment this out and use 
         ! the checkpointed variables xxx_fwd. 
!        do idow2 = 1, idow-1
!          gasqd = gamaan*gamaan
!          wsqd = wh2o*wh2o
!          kw2 = kan*wsqd/gasqd
!          aa = 1.-kw2
!          bb = twoso4+kw2*(tno3+tnh4-twoso4)
!          cc = -(kw2*tno3*(tnh4-twoso4))
!          disc = bb*bb-4.*aa*cc
!          if (aa .ne. 0.d0) then
!            dd = sqrt(disc)
!            xxq = -(0.5d0*(bb+sign(1.d0,bb)*dd))
!            rr1 = xxq/aa
!            rr2 = cc/xxq
!            if (rr1*rr2 .lt. 0.d0) then
!              xno3 = max(rr1,rr2)
!            else
!              xno3 = min(rr1,rr2)
!            endif
!          else
!            xno3 = -(cc/bb)
!          endif
!          xno3 = min(xno3,tno3)
!          call awater( irh,tso4,ynh4,xno3,ah2o )
!          wh2o = 0.001*ah2o
!          man = xno3/wh2o
!          mas = tso4/wh2o
!          mnh4 = 2.*mas+man
!          ynh4 = mnh4*wh2o
!          cat(1) = 0.
!          cat(2) = mnh4
!          an(1) = mas
!          an(2) = man
!          an(3) = 0.
!          call actcof( cat,an,gams,molnu,phibar )
!          gamaan = gams(2,2)
!          gamold = gamaan
!        end do
        IF ( idow .gt. 1 ) THEN
           gamaan =  gamaan_fwd(I,J,L,idow-1)
           wh2o   =  wh2o_fwd(I,J,L,idow-1)
           ynh4   =  ynh4_fwd(I,J,L,idow-1)
        ENDIF
        gamold = gamaan
        ynh4i = ynh4
        wh2oi = wh2o
        gamaani = gamaan
        gasqd = gamaan*gamaan
        wsqd = wh2o*wh2o
        kw2 = kan*wsqd/gasqd
        aa = 1.-kw2
        bb = twoso4+kw2*(tno3+tnh4-twoso4)
        cc = -(kw2*tno3*(tnh4-twoso4))
        disc = bb*bb-4.*aa*cc
        if (aa .ne. 0.d0) then
          dd = sqrt(disc)
          xxq = -(0.5d0*(bb+sign(1.d0,bb)*dd))
          rr1 = xxq/aa
          rr2 = cc/xxq
          if (rr1*rr2 .lt. 0.d0) then
            xno3 = max(rr1,rr2)
          else
            xno3 = min(rr1,rr2)
          endif
        else
          xno3 = -(cc/bb)
        endif
        xno3 = min(xno3,tno3)
        call awater( irh,tso4,ynh4,xno3,ah2o )
        wh2o = 0.001*ah2o
        man = xno3/wh2o
        mas = tso4/wh2o
        mnh4 = 2.*mas+man
        ynh4 = mnh4*wh2o
        cat(1) = 0.
        cat(2) = mnh4
        an(1) = mas
        an(2) = man
        an(3) = 0.
        call actcof( cat,an,gams,molnu,phibar )
        gamaan = gams(2,2)
        eror = abs(gamold-gamaan)/gamold
        if (eror .le. toler1) then
          ano3 = xno3*mwno3
          adwh2o = adwh2o+1000*adah2o
          adah2o = 0.
          adtnh4 = adtnh4+adgnh3*mwnh3*(0.5-sign(0.5,floor-(tnh4-ynh4)))
          adynh4 = adynh4-adgnh3*mwnh3*(0.5-sign(0.5,floor-(tnh4-ynh4)))
          adgnh3 = 0.
          adano3 = adano3-adgno3*(0.5-sign(0.5,floor-(tmasshno3-ano3)))
          adtmasshno3 = adtmasshno3+adgno3*(0.5-sign(0.5,floor-
     $(tmasshno3-ano3)))
          adgno3 = 0.
          adynh4 = adynh4+adanh4*mwnh4
          adanh4 = 0.
          adxno3 = adxno3+adano3*mwno3
          adano3 = 0.
          adtso4 = adtso4+adaso4*mwso4
          adaso4 = 0.
        endif
        adgamaan = adgamaan+adgamold
        adgamold = 0.
        aderori = aderor/gamold
        adgamold = adgamold-aderor*(abs(gamold-gamaan)/(gamold*gamold))
        adgamaan = adgamaan-aderori*sign(1.,gamold-gamaan)
        adgamold = adgamold+aderori*sign(1.,gamold-gamaan)
        aderor = 0.
        adgams(2,2) = adgams(2,2)+adgamaan
        adgamaan = 0.
        call adactcof( cat,an,adcat,adan,adgams )
        adan(3) = 0.
        adman = adman+adan(2)
        adan(2) = 0.
        admas = admas+adan(1)
        adan(1) = 0.
        admnh4 = admnh4+adcat(2)
        adcat(2) = 0.
        adcat(1) = 0.
        admnh4 = admnh4+adynh4*wh2o
        adwh2o = adwh2o+adynh4*mnh4
        adynh4 = 0.
        adman = adman+admnh4
        admas = admas+2*admnh4
        admnh4 = 0.
        adtso4 = adtso4+admas/wh2o
        adwh2o = adwh2o-admas*(tso4/(wh2o*wh2o))
        admas = 0.
        adwh2o = adwh2o-adman*(xno3/(wh2o*wh2o))
        adxno3 = adxno3+adman/wh2o
        adman = 0.
        adah2o = adah2o+0.001*adwh2o
        adwh2o = 0.
        ynh4 = ynh4i
        call adawater( irh,tso4,ynh4,xno3,adtso4,adynh4,adxno3,adah2o )
        if (aa .ne. 0.d0) then
          dd = sqrt(disc)
          xxq = -(0.5d0*(bb+sign(1.d0,bb)*dd))
          rr1 = xxq/aa
          rr2 = cc/xxq
          if (rr1*rr2 .lt. 0.d0) then
            xno3 = max(rr1,rr2)
          else
            xno3 = min(rr1,rr2)
          endif
        else
          xno3 = -(cc/bb)
        endif
        adtno3 = adtno3+adxno3*(0.5-sign(0.5,tno3-xno3))
        adxno3 = adxno3*(0.5+sign(0.5,tno3-xno3))
        if (aa .ne. 0.d0) then
          if (rr1*rr2 .lt. 0.d0) then
            adrr1 = adrr1+adxno3*(0.5+sign(0.5,rr1-rr2))
            adrr2 = adrr2+adxno3*(0.5-sign(0.5,rr1-rr2))
            adxno3 = 0.
          else
            adrr1 = adrr1+adxno3*(0.5+sign(0.5,rr2-rr1))
            adrr2 = adrr2+adxno3*(0.5-sign(0.5,rr2-rr1))
            adxno3 = 0.
          endif
          adcc = adcc+adrr2/xxq
          adxxq = adxxq-adrr2*(cc/(xxq*xxq))
          adrr2 = 0.
          adaa = adaa-adrr1*(xxq/(aa*aa))
          adxxq = adxxq+adrr1/aa
          adrr1 = 0.
          adbb = adbb-0.5d0*adxxq
          addd = addd-0.5d0*adxxq*sign(1.d0,bb)
          adxxq = 0.
          addisc = addisc+addd*(1./(2.*sqrt(disc)))
          addd = 0.
        else
          adbb = adbb+adxno3*(cc/(bb*bb))
          adcc = adcc-adxno3/bb
          adxno3 = 0.
        endif
        adaa = adaa-4*addisc*cc
        adbb = adbb+2*addisc*bb
        adcc = adcc-4*addisc*aa
        addisc = 0.
        adkw2 = adkw2-adcc*tno3*(tnh4-twoso4)
        adtnh4 = adtnh4-adcc*kw2*tno3
        adtno3 = adtno3-adcc*kw2*(tnh4-twoso4)
        adtwoso4 = adtwoso4+adcc*kw2*tno3
        adcc = 0.
        adkw2 = adkw2+adbb*(tno3+tnh4-twoso4)
        adtnh4 = adtnh4+adbb*kw2
        adtno3 = adtno3+adbb*kw2
        adtwoso4 = adtwoso4+adbb*(1-kw2)
        adbb = 0.
        adkw2 = adkw2-adaa
        adaa = 0.
        adgasqd = adgasqd-adkw2*(kan*wsqd/(gasqd*gasqd))
        adwsqd = adwsqd+adkw2*(kan/gasqd)
        adkw2 = 0.
        wh2o = wh2oi
        adwh2o = adwh2o+2*adwsqd*wh2o
        adwsqd = 0.
        gamaan = gamaani
        adgamaan = adgamaan+2*adgasqd*gamaan
        adgasqd = 0.
      end do
      adtwoso4 = adtwoso4+adynh4
      adynh4 = 0.
      adah2o = adah2o+0.001d0*adwh2o
      adwh2o = 0.
      ynh4 = twoso4
      call adawater( irh,tso4,ynh4,tno3,adtso4,adynh4,adtno3,adah2o )
      adtwoso4 = adtwoso4+adynh4
      adynh4 = 0.
      adtso4 = adtso4+2*adtwoso4
      adtwoso4 = 0.
      ano3 = in(4)
      adano3 = adano3+adtmasshno3*(0.5-sign(0.5d0,0.d0-(gno3+ano3)))
      adgno3 = adgno3+adtmasshno3*(0.5-sign(0.5d0,0.d0-(gno3+ano3)))
      adtmasshno3 = 0.
      anh4 = in(5)
      adanh4 = adanh4+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh4)
      adgnh3 = adgnh3+adtnh4*((0.5-sign(0.5d0,0.d0-(gnh3/mwnh3+anh4/
     $mwnh4)))/mwnh3)
      adtnh4 = 0.
      adano3 = adano3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwno3)
      adgno3 = adgno3+adtno3*((0.5-sign(0.5d0,0.d0-(ano3/mwno3+gno3/
     $mwhno3)))/mwhno3)
       adtno3 = 0.
      adso4 = adso4+adtso4*((0.5-sign(0.5,floor-so4/mwso4))/mwso4)
      adtso4 = 0.
      adin(5) = adin(5)+adanh4
      adanh4 = 0.
      adin(4) = adin(4)+adano3
      adano3 = 0.
      adin(3) = adin(3)+adgnh3
      adgnh3 = 0.
      adin(2) = adin(2)+adgno3
      adgno3 = 0.
      adin(1) = adin(1)+adso4
      adso4 = 0.


      end SUBROUTINE ADRPMARES_6_D5

!------------------------------------------------------------------------------

      ! End of module
      END MODULE RPMARES_ADJ_MOD   
