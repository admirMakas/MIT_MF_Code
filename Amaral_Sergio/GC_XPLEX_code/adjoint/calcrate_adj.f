! $Id: calcrate_adj.f,v 1.1 2010/04/01 07:09:43 daven Exp $
!      SUBROUTINE CALCRATE( SUNCOS )
       SUBROUTINE CALCRATE_ADJ(RRATE_ADJ,IX,IY,IZ)

!
!******************************************************************************
!  Subroutine CALCRATE_ADJ basically just transfers adjoints of emissions and
!  deposition rates from the RRATE_ADJ array to more specific arrays.  This 
!  is only for species whose emission and/or deposition is handled within the
!  fullchemistry mechanims (such as NOx, but not SOx). (dkh, 06/05/06)  
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) RRATE_ADJ (REAK*8(:)) : Adjoint of emission and deposition rates
!  (2-4) IX, IY, IZ (INTEGER) : 3-D array location             [unit]
!     
!  Module variable as Input:
!  ============================================================================
!     
!  Module variable as Output:
!  ============================================================================
!  (1 ) ADJ_REMIS  (REAK*8(:,:))   : Adjoint of emission rates
!  (2 ) ADJ_DEPSAV (REAK*8(:,:,:)) : Adjoint of deposition rates
!  (3 ) ADJ_TAREA  (REAK*8(:,:))   : Adjoint of areosol area
!
!  NOTES:
!  (1 ) Updated to GCv8 (dkh, 03/30/10) 
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE COMODE_MOD,     ONLY : ABSHUM, AIRDENS, ERADIUS, T3, TAREA
      USE COMODE_MOD,     ONLY : JLOP
      USE DRYDEP_MOD,     ONLY : NUMDEP 
      USE ERROR_MOD,      ONLY : ERROR_STOP, GEOS_CHEM_STOP
      USE PBL_MIX_MOD,    ONLY : GET_FRAC_UNDER_PBLTOP
      
      USE ADJ_ARRAYS_MOD, ONLY : DEPSAV_ADJ, REMIS_ADJ
      USE GCKPP_ADJ_GLOBAL, ONLY : NCOEFF
 
      USE MYTYPE 
      USE COMPLEXIFY 
      IMPLICIT NONE 

#     include "CMN_SIZE"    ! Size parameters
#     include "comode.h"    ! NTDEP, NEMIS, NTEMIS, NCSURBAN


      ! Arguments
      TYPE (XPLEX),INTENT(IN) :: RRATE_ADJ(NMTRATE)
      INTEGER         :: IX, IY, IZ
 
      ! Local variables 
      INTEGER :: NK, I, NN, N
      INTEGER :: JJLOOP
      TYPE (XPLEX)  :: XAREA,XRADIUS,XSQM
      TYPE (XPLEX)  :: XSTKCF,XDENA,XSTK
      TYPE (XPLEX)  :: ARSL1K, DFKG

      !TYPE (XPLEX)  :: ADJ_ARSL1K, ADJ_XAREA
      !TYPE (XPLEX)  :: ADJ_XRADIUS
      TYPE (XPLEX)  :: ARSL1K_ADJ, XAREA_ADJ
      TYPE (XPLEX)  :: XRADIUS_ADJ
      !=================================================================
      !CALCRATE_ADJ begins here!
      !=================================================================
C
C *********************************************************************
C ******         ADJOINT OF SET DRY DEPOSITION RATES             ******
C *********************************************************************
C
      DO I = 1,NUMDEP 
         NK = NTDEP(I)
         IF (NK.NE.0) THEN
            ! We don't loop over SMVG blocks for adjoint
            !DO KLOOP = 1,KTLOOP

               ! Pass JJLOOP, IX, IY and IZ as arguments instead (dkh, 06/04/06)  
               ! 1-D grid box index (accounts for reordering)
               !JLOOP = LREORDER(KLOOP+JLOOPLO)

               ! 3-D grid box index
               !IX    = IXSAVE(JLOOP)
               !IY    = IYSAVE(JLOOP)
               !IZ    = IZSAVE(JLOOP)

               ! Adjoint of deposition frequency 
               ! PBLFRAC is the fraction of grid box (I,J,L) below the PBL top
               ! fwd code:
               !RRATE(KLOOP,NK) = DEPSAV(IX,IY,I) *      
               !             GET_FRAC_UNDER_PBLTOP( IX, IY, IZ )

               DEPSAV_ADJ(IX,IY,I) 
     &            = RRATE_ADJ(NK)  * 
     &                           GET_FRAC_UNDER_PBLTOP( IX, IY, IZ )


            !ENDDO
         ENDIF
      ENDDO

C *********************************************************************
C ******            ADJOINT OF SET EMISSION RATES                ******
C *********************************************************************
C

      NCS = 1 
      DO I = 1,NEMIS(NCS)
C get tracer number corresponding to emission species I
         NN = IDEMS(I)
         IF (NN.NE.0) THEN
C find reaction number for emission of tracer NN
            NK = NTEMIS(NN,NCS)
            IF (NK.NE.0) THEN
               ! We don't loop over SMVG blocks for adjoint
               !DO KLOOP = 1,KTLOOP
               
                  ! Pass JJLOOP as an argument for adjoint
                  !JLOOP = LREORDER(KLOOP+JLOOPLO)

                  ! fwd code:
                  ! RRATE(KLOOP,NK) = REMIS(JLOOP,I)
                  ! At this point, all the adjoint routine has to do is pass the 
                  ! values from RRATE_ADJ to REMIS_ADJ. 
                  JJLOOP              = JLOP(IX,IY,IZ) 
                  REMIS_ADJ(JJLOOP,I) = RRATE_ADJ(NK)

               !ENDDO
            ENDIF
         ENDIF
      ENDDO

! aerosol het chem adjoint, need to update
!
!      NCS = NCSURBAN
!
!        ! Set HETCHEM = T to perform het chem on aerosols
!        !HETCHEM = .TRUE.
!
!        !IF ( HETCHEM ) THEN
!
!           ! Initialize TAREA_ADJ
!           TAREA_ADJ(JJLOOP,:) = 0d0
!
!           !===========================================================
!           ! Perform heterogeneous chemistry on sulfate aerosol
!           ! plus each of the NDUST dust size bins from FAST-J
!           !===========================================================
!           XDENA   = AIRDENS(JJLOOP)
!           XSTK    = SQRT(T3(JJLOOP))
!
!           DO I       = 1, NNADDK(NCS)
!              NK      = NKSPECK(I,NCS)
!              XSQM    = SQRT(ARR(NK,NCS))
!
!              ARSL1K_ADJ          = 0d0
!
!              ! Loop over sulfate and other aerosols
!              ! SKIPP DUST for now
!              !DO N = 1, NDUST + NAER
!              !DO N = NDUST+1, NDUST + NAER
!              ! Now include carbon aerosol (dkh, 06/03/08) 
!              DO N = NDUST+1, NDUST + 3
!
!                 ! Adjoint of ARSL1K
!                 ! fwd code: 
!                 !RRATE(KLOOP,NK) = RRATE(KLOOP,NK) + ARSL1K  
!                 ARSL1K_ADJ = RRATE_ADJ(JJLOOP,NK)
!
!                 ! Recalculate XSTKCF, XSQM, XRADIUS, XAREA
!
!                 ! Surface area of aerosol [cm2 aerosol/cm3 air]
!                 XAREA = TAREA(JJLOOP,N)
!
!                 ! Test if N2O5 hydrolysis rxn
!                 IF ( NK == NKN2O5 ) THEN
!
!                    ! Get GAMMA for N2O5 hydrolysis, which is
!                    ! a function of aerosol type, temp, and RH
!                    XSTKCF = N2O5( N, T3(JJLOOP), ABSHUM(JJLOOP) )
!
!                 ELSE
!
!                    ! Get GAMMA for species other than N2O5
!                    XSTKCF = BRR(NK,NCS)
!
!                 ENDIF
!
!                 ! Radius for dust size bin N
!                 XRADIUS = ERADIUS(JJLOOP,N)
!
!
!                 ! ARSL1K begins here!
!                 !=================================================================
!                 IF ( XAREA < 0d0 .or. XRADIUS < 1d-30 ) THEN
!
!                    ! fwd code
!                    !ARSL1K = 1.D-3
!                    ! Adjoint of this is do nothing
!                    XAREA_ADJ   = 0d0 
!                    XRADIUS_ADJ = 0d0 
!
!                 ELSE
!
!                    ! Recalculate DFKG 
!                    ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
!                    DFKG  = 9.45D17/XDENA * XSTK 
!     &                    * SQRT(3.472D-2 + 1.D0/(XSQM*XSQM))
!      
!                    ! Calcualte adjoint of AREA from ARSL1K_ADJ 
!                    ! fwd code 
!                    !ARSL1K = AREA / ( RADIUS/DFKG + 2.749064E-4*SQM/(STKCF*STK) )
!                    XAREA_ADJ = ARSL1K_ADJ / ( XRADIUS/DFKG 
!     &                        + 2.749064E-4*XSQM/(XSTKCF*XSTK) )
!
!                    ! Calculate adjoint of RADIUS from ARSL1K_ADJ
!                    XRADIUS_ADJ = ARSL1K_ADJ * ( - XAREA / DFKG ) 
!     &                          * ( XRADIUS/DFKG
!     &                          + 2.749064E-4*XSQM/(XSTKCF*XSTK) ) ** -2
!
!
!                 ENDIF
! 
!                 ! Surface area of aerosol [cm2 aerosol/cm3 air]
!                 ! fwd code:
!                 !XAREA = TAREA(JLOOP,N)
!                 TAREA_ADJ(JJLOOP,N)  = TAREA_ADJ(JJLOOP,N) + XAREA_ADJ
!
!                 ! fwd code:
!                 !XRADIUS = ERADIUS(JLOOP,N) 
!                 ERADIUS_ADJ(JJLOOP,N) = ERADIUS_ADJ(JJLOOP,N) 
!     &                                + XRADIUS_ADJ
!
!              ENDDO
!
!              ! Reset, not needed, but to be safe...
!              !RRATE_ADJ(JJLOOP,NK) = 0d0
!
!            ENDDO
!  
!         !ENDIF 
!      !ENDDO 
C
      RETURN

C
C *********************************************************************
C                       INTERNAL SUBROUTINES 
C *********************************************************************
C
      CONTAINS

      FUNCTION N2O5( AEROTYPE, TEMP, RH ) RESULT( GAMMA )

      !=================================================================
      ! Internal function N2O5 computes the GAMMA sticking factor
      ! for N2O5 hydrolysis. (mje, bmy, 8/7/030
      ! 
      ! Arguments as Input:
      ! ----------------------------------------------------------------
      ! (1 ) AEROTYPE (INTEGER) : # denoting aerosol type (cf FAST-J)
      ! (2 ) TEMP     (TYPE (XPLEX) ) : Temperature [K]
      ! (3 ) RH       (TYPE (XPLEX) ) : Relative Humidity [fraction]
      !
      ! NOTES:
      !=================================================================

      ! Arguments
      INTEGER, INTENT(IN) :: AEROTYPE
      TYPE (XPLEX),  INTENT(IN) :: TEMP, RH

      ! Local variables
      TYPE (XPLEX)              :: RH_P, FACT, TTEMP

      ! Function return value
      TYPE (XPLEX)              :: GAMMA

      !=================================================================
      ! N2O5 begins here!
      !=================================================================

      ! Convert RH to % (max = 100%)
      RH_P  = MIN( RH * 100d0, 100d0 )

      ! Default value
      GAMMA = 0.01d0

      ! Special handling for various aerosols
      SELECT CASE ( AEROTYPE )

         !----------------
         ! Dust 
         !----------------
         CASE ( 1, 2, 3, 4, 5, 6, 7 )

            ! Based on unpublished Crowley work
            GAMMA = 0.01d0

         !----------------
         ! Sulfate
         !----------------
         CASE ( 8 )

            !===========================================================
            ! RH dependence from Kane et al., Heterogenous uptake of 
            ! gaseous N2O5 by (NH4)2SO4, NH4HSO4 and H2SO4 aerosols
            ! J. Phys. Chem. A , 2001, 105, 6465-6470 
            !===========================================================
            GAMMA = 2.79d-4 + RH_P*(  1.30d-4 +
     &                        RH_P*( -3.43d-6 +
     &                        RH_P*(  7.52d-8 ) ) )
 
            !===========================================================
            ! Temperature dependence factor (Cox et al, Cambridge UK) 
            ! is of the form:
            !
            !          10^( LOG10( G294 ) - 0.04 * ( TTEMP - 294 ) )
            ! FACT = -------------------------------------------------
            !                     10^( LOG10( G294 ) )
            !
            ! Where G294 = 1e-2 and TTEMP is MAX( TEMP, 282 ).
            ! 
            ! For computational speed, replace LOG10( 1e-2 ) with -2
            ! and replace 10^( LOG10( G294 ) ) with G294 
            !===========================================================
            TTEMP = MAX( TEMP, 282d0 )
            FACT  = 10.d0**( -2d0 - 4d-2*( TTEMP - 294.d0 ) ) / 1d-2

            ! Apply temperature dependence
            GAMMA = GAMMA * FACT

         !----------------
         ! Black Carbon
         !----------------
         CASE ( 9 )

             ! From IUPAC
             GAMMA = 0.005d0

         !----------------
         ! Organic Carbon
         !----------------           
         CASE ( 10 )

            !===========================================================
            ! Based on Thornton, Braban and Abbatt, 2003
            ! N2O5 hydrolysis on sub-micron organic aerosol: the effect
            ! of relative humidity, particle phase and particle size
            !===========================================================
            IF ( RH_P >= 57d0 ) THEN
               GAMMA = 0.03d0
            ELSE
               GAMMA = RH_P * 5.2d-4
            ENDIF

         !----------------
         ! Sea salt
         ! accum & coarse
         !----------------
         CASE ( 11, 12 )

            ! Based on IUPAC recomendation
            IF ( RH_P >= 62 ) THEN
               GAMMA = 0.03d0
            ELSE
               GAMMA = 0.005d0
            ENDIF

         !----------------         
         ! Default
         !----------------
         CASE DEFAULT
            WRITE (6,*) 'Not a suitable aerosol surface '
            WRITE (6,*) 'for N2O5 hydrolysis'
            WRITE (6,*) 'AEROSOL TYPE =',AEROTYPE
            CALL GEOS_CHEM_STOP

      END SELECT

      ! Return to CALCRATE
      END FUNCTION N2O5 

      ! Return to calling program
      END SUBROUTINE CALCRATE_ADJ
                                 
