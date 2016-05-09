!$Id: linoz_adj_mod.f,v 1.6 2012/05/08 02:18:25 nicolas Exp $
      MODULE LINOZ_ADJ_MOD

C Revision 2.10  2000/03/23 20:41:45  pjc
C Initial version adapted heavily from McLinden's original file.

c
c----------------------------------------------------------------------- 

      USE MYTYPE
      USE COMPLEXIFY
      CONTAINS

      !============================================================

      subroutine do_linoz_adj

      USE TIME_MOD

#     include "CMN_SIZE"

      ! Local variables
      ! Chem time step in seconds for linoz   (dbj,bdf 6/24/03)
      TYPE (XPLEX)  :: NSCHEM

      NSCHEM = GET_TS_CHEM()*60D0 ! Linoz needs time step in seconds
      CALL LINOZ_CHEM3_ADJ(NSCHEM)

      end subroutine do_linoz_adj   

!-------------------------------------------------------------------
      SUBROUTINE LINOZ_CHEM3_ADJ( DTCHEM )

!
!***************************************************************
! Subroutine LINOZ_CHEM3_ADJ is the adjont of LINOZ_CHEM3, 
! manually derived to account for strat flux adjoints.
! 
! This replaces an older version of this routine that was 
! based on TAMC code.  ( hml, dkh, 02/20/12, adj32_025)
!
!***************************************************************
      USE TIME_MOD,        ONLY : GET_NHMS
      USE TIME_MOD,        ONLY : GET_NYMD
      USE DAO_MOD,         ONLY : AD
      USE DAO_MOD,         ONLY : T
      USE ERROR_MOD,       ONLY : ERROR_STOP
      USE TRACER_MOD,      ONLY : TCVV      
      USE TRACER_MOD,      ONLY : STT_TMP
      USE TRACER_MOD,      ONLY : ITS_A_TAGOX_SIM
      USE TRACERID_MOD,    ONLY : IDTOX
      USE GRID_MOD,        ONLY : GET_AREA_CM2   
      USE TROPOPAUSE_MOD,  ONLY : GET_TPAUSE_LEVEL
      USE TROPOPAUSE_MOD,  ONLY : GET_MAX_TPAUSE_LEVEL
      USE PRESSURE_MOD,    ONLY : GET_PEDGE
      USE PRESSURE_MOD,    ONLY : GET_PCENTER
      USE CHECKPOINT_MOD,  ONLY : READ_UPBDFLX_CHKFILE
      USE LOGICAL_ADJ_MOD, ONLY : LADJ_STRAT
      USE ADJ_ARRAYS_MOD,  ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : PROD_SF,     LOSS_SF
      USE ADJ_ARRAYS_MOD,  ONLY : PROD_SF_ADJ, LOSS_SF_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD, LFD
      USE ADJ_ARRAYS_MOD,  ONLY : Ox_p, Ox_l

  
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "CMN"
!#     include "../new/linoz.com"

      TYPE (XPLEX), INTENT(IN) :: DTCHEM    ! Time step [seconds]


C==============================================
C define arguments
C==============================================

      ! hml: add for strat prod & loss sense
      INTEGER :: IMX,     JM,         LM
      INTEGER :: I,       J,          L
      INTEGER :: LBOT,    L_OVERWRLD
      INTEGER :: NTRACER, NUM_TRACER, LPOS,   ITRC
      INTEGER :: NHMS,    NYMD

      TYPE (XPLEX)  :: CLIMO3,  CLIMPML,    PMLTOT
      TYPE (XPLEX)  :: DCO3,    DERO3,      DERTMP
      TYPE (XPLEX)  :: DERCO3,  DMASS,      DTMP
      TYPE (XPLEX)  :: SSO3

      TYPE (XPLEX)  :: TAU
      TYPE (XPLEX)  :: P,        k,       M0
      TYPE (XPLEX)  :: P_ADJ,    k_ADJ,   M0_ADJ
      TYPE (XPLEX)  :: LOSS_ADJ, PROD_ADJ
      TYPE (XPLEX)  :: PROD,     LOSS
      TYPE (XPLEX)  :: PROD_0,   LOSS_0

      ! Arrays
      TYPE (XPLEX)  :: DCOLO3(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)  :: COLO3(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)  :: OUT_DATA(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)  :: TLSTT(JJPAR,LLPAR,7)

      TYPE (XPLEX)  :: STT_ADJ_IN(IIPAR,JJPAR,LLPAR)

      ! Assign values for local IMX and JM   (dbj 6/24/03) 
      IMX    = IIPAR
      JM     = JJPAR
      LM     = LLPAR   ! dbj

      L_OVERWRLD = GET_MAX_TPAUSE_LEVEL()

      NTRACER = IDTOX

      ! READING STT and TLSTT IN REVERSE MODE       
      NHMS     = GET_NHMS()
      NYMD     = GET_NYMD()
      CALL READ_UPBDFLX_CHKFILE( NYMD, NHMS )

      DO L = 1,LLPAR
      DO J = 1,JJPAR
      DO I = 1,7
         TLSTT(J,L,I) = STT_TMP(I,J,L,1)
      ENDDO
      ENDDO
      ENDDO

      ! don't overwrite (hml, 10/15/11)
      !STT(:,:,:,NTRACER) = STT_TMP(:,:,:,2)

      WRITE(6,*) '-----------------------------------------------------'
      write(6,*) '     doing adjoint linoz stratospheric chemistry     '
      WRITE(6,*) '-----------------------------------------------------'

      STT_ADJ_IN(:,:,:) = STT_ADJ(:,:,:,IDTOX)

      !OUT_DATA = 0d0

      ! Initialize arrays (hml, 10/17/11)
      LOSS = 0d0
      PROD = 0d0

      ! NEW: 
      ! Now use this format at all times (dkh, 04/20/12) 
      !! For strat P & L optimization (hml, 10/03/11)
      DO J = 1, JM
      DO I = 1, IMX
         LBOT = GET_TPAUSE_LEVEL(I,J)+1
         LPOS = 1
   
         ! To set LFD properly (hml, 10/12/11)
         IF ( I == IFD.and.J == JFD) THEN 
            print *, 'LBOT = ', LBOT
         ENDIF
   
         DO WHILE (GET_PEDGE(I,J,LPOS+1) .GE. 0.3D0)
            LPOS = LPOS +1
         ENDDO
         LPOS = LPOS-1
             
         !---------------------------------------------------------
         ! dbj: for now, set tagged stratospheric tracer to total 
         ! O3 in the overworld to avoid issues with spin ups
         !---------------------------------------------------------
         IF ( ITS_A_TAGOX_SIM() ) THEN
           STT_TMP(I,J,(L_OVERWRLD+1):LLPAR,NTRACER) =
     &              STT_TMP(I,J,(L_OVERWRLD+1):LLPAR,1)
         ENDIF

! With this format we need to start at LLPAR so that COLO3 and DCOLO3 are correct. 
!!!            ! If we just loop from LPOS, rather than LLPAR, then we only deal with
!!!            ! levels for which PEDGE > 0.3d0
         DO L = LM,LBOT,-1        

            IF (STT_TMP(I,J,L,2) .LE. 0.D0) CYCLE

            !---------------------------------------
            ! GET RATES - assigning PROD and LOSS
            !---------------------------------------

            !------------------------------------------------------
            ! Recalculate forward model values to get rates
            !------------------------------------------------------

            ! bdf stt is in v/v, make conversion to DU
            IF ( L .EQ. LM) THEN !top model layer
               ! Use checkpointed value
               !DCOLO3(I,J,L) = (STT(I,J,L,NTRACER)*AD(I,J,L)/
               DCOLO3(I,J,L) = (STT_TMP(I,J,L,2)*AD(I,J,L)/
     &            TCVV(NTRACER))/ GET_AREA_CM2(J) *
     &            6.022d23/(28.97/TCVV(NTRACER)*1d-3)/ 2.687d16
               COLO3(I,J,L) = DCOLO3(I,J,L)*0.5
            ELSE
               ! Use checkpointed value
               !DCOLO3(I,J,L) = (STT(I,J,L,NTRACER)*AD(I,J,L)/
               DCOLO3(I,J,L) = (STT_TMP(I,J,L,2)*AD(I,J,L)/
     &              TCVV(NTRACER))/ GET_AREA_CM2(J) *
     &              6.022d23/(28.97/TCVV(NTRACER)*1d-3)/ 2.687d16
               COLO3(I,J,L) = COLO3(I,J,L+1) +
     &               (DCOLO3(I,J,L)+DCOLO3(I,J,L+1))*0.5
            ENDIF

            ! ++++++ climatological P-L:   ++++++          
            CLIMPML = TLSTT(J,L,4)      ! Climatological P-L = (P-L)^o
     
            ! ++++++ local ozone feedback: ++++++ 
            DERO3 = TLSTT(J,L,5)      ! Derivative w.r.t. O3. dero3=-1/(time constant)
            IF ( DERO3 .EQ. 0 ) CYCLE ! Skip Linoz if lifetime is infinite.
            CLIMO3 = TLSTT(J,L,1)     ! Climatological O3 = f^o
            DERCO3 = TLSTT(J,L,7)     ! Derivative w.r.t. Column O3
            DCO3   = (COLO3(I,J,L) - TLSTT(J,L,3)) ! deviation from o3 climatology.
    
            ! ++++++ temperature feedback: ++++++ 
            DERTMP = TLSTT(J,L,6)              ! Derivative w.r.t. Temperature
            DTMP   = (T(I,J,L) - TLSTT(J,L,2)) ! Deviation in Temperature from climatology.
   
            ! ++++++ calculate steady-state ozone: ++++++ 
            SSO3 = CLIMO3 
     &           - (CLIMPML + DTMP*DERTMP + DCO3*DERCO3) /DERO3
      
            ! debug: recalculated DMASS just to check with fwd model
            ! Use checkpointed value
            !DMASS = (SSO3 - STT(I,J,L,NTRACER))
            DMASS = (SSO3 - STT_TMP(I,J,L,2))
     &            * (1.0 - exp(DERO3*DTCHEM))
      
            ! note: there is a factor of TC / AD * AD / TC that cancels 
            ! out in definition of PROD_0
            PROD_0 = - (SSO3 * DERO3) 
            LOSS_0 = - DERO3

            !---------------------------------------
            ! END of GET RATES 
            !---------------------------------------

            IF (GET_PEDGE(I,J,L) .LE. 0.3D0) THEN

               ! fwd:
               !STT(I,J,L,NTRACER) = ( GET_PCENTER(I,J,L)
               !                   / GET_PCENTER(I,J,LPOS-1) )
               !                   * STT(I,J,LPOS-1,NTRACER)
               ! adj:
               STT_ADJ(I,J,LPOS-1,NTRACER) 
     &            = STT_ADJ(I,J,LPOS-1,NTRACER)
     &            + ( GET_PCENTER(I,J,L) / GET_PCENTER(I,J,LPOS-1) )
     &            * STT_ADJ(I,J,L,NTRACER) 
               

            !otherwise just take the adjoint of the low-pressure decay
            ! and the prod / loss scaling factors have no effect
            ELSE
    
               !! Scaled prod & loss rate   
               IF ( LADJ_STRAT ) THEN 
                  PROD = PROD_0 * PROD_SF(I,J,1,Ox_p)
                  LOSS = LOSS_0 * LOSS_SF(I,J,1,Ox_l)
               ENDIF 

               k = LOSS                             ! loss freq [s-1]
               P = PROD * AD(I,J,L) / TCVV(NTRACER) ! production term [kg s-1]
               ! Use checkpointed value
               ! Put ozone back to kg (hml, 11/06/11) 
               M0 = STT_TMP(I,J,L,2) 
     &            * AD(I,J,L) / TCVV(NTRACER)! initial mass [kg]

               ! No prod or loss at all
               if ( k .eq. 0d0 .and. P .eq. 0d0 ) cycle

               ! fwd code:
               !STT(I,J,L,N) = M0 * exp(-k*t) + (P/k)*(1d0-exp(-k*t))
               ! adj code: note use the input value of STT_ADJ and 
               ! convert units of STT_ADJ to pre LINOZ_ADJ units
               M0_ADJ = STT_ADJ_IN(I,J,L) * TCVV(NTRACER) / AD(I,J,L)
     &                * exp(-k*DTCHEM)
               P_ADJ  = STT_ADJ_IN(I,J,L) * TCVV(NTRACER) / AD(I,J,L)
     &                * (1d0 - exp(-k*DTCHEM))/k
               k_ADJ  = STT_ADJ_IN(I,J,L) * TCVV(NTRACER) / AD(I,J,L) 
     &                * ( -p/(k**2) + p/(k**2)*exp(-k*DTCHEM)
     &                + (p*DTCHEM/k)*exp(-k*DTCHEM) 
     &                - DTCHEM * exp(-k*DTCHEM) * M0 )


               ! fwd code:
               !k = LOSS(I,J,L,N)                       ! loss freq [s-1]
               !P = PROD(I,J,L,N) * AD(I,J,L) / TCVV(N) ! production term [kg s-1]
               !M0 = STT(I,J,L,N)                       ! initial mass [kg]
               ! adj code:
               LOSS_ADJ          = K_ADJ
               PROD_ADJ          = P_ADJ * AD(I,J,L) / TCVV(NTRACER)

               !!! Now calculate the update to STT_ADJ here.
               STT_ADJ (I,J,L,NTRACER) = M0_ADJ / TCVV(NTRACER) 
     &                                    * AD(I,J,L)

               !------------------------------------------------------
               ! adjoint with respect to PROD and LOSS scaling factors
               !------------------------------------------------------
               IF ( LADJ_STRAT ) THEN 
                  ! fwd code:
                  !PROD(I,J,L,N) = PROD_0(I,J,L,N) * PROD_SF(I,J,1,N)
                  !LOSS(I,J,L,N) = LOSS_0(I,J,L,N) * LOSS_SF(I,J,1,N)
                  ! adj code:
                  PROD_SF_ADJ(I,J,1,Ox_p) = PROD_SF_ADJ(I,J,1,Ox_p)
     &                                 + PROD_0 * PROD_ADJ
                  LOSS_SF_ADJ(I,J,1,Ox_l) = LOSS_SF_ADJ(I,J,1,Ox_l)
     &                                 + LOSS_0 * LOSS_ADJ
               ENDIF 

            ENDIF  ! PEDGE

         ENDDO       ! loop over L

      ENDDO          ! loop over I
      ENDDO          ! loop pver J


      !write our calculated column o3 maximum
      !write(6,*) 'max of columns= ',maxval(out_data)

!!$OMP END PARALLEL DO 

      END SUBROUTINE LINOZ_CHEM3_ADJ 
!------------------------------------------------------------------------------

      ! End of module
      END MODULE LINOZ_ADJ_MOD
