! $Id: partition.f,v 1.2 2012/03/01 22:00:27 daven Exp $
      SUBROUTINE PARTITION( NTRACER, STT, XNUMOL ) 
!
!******************************************************************************
!  Subroutine PARTITION separates GEOS-CHEM tracers into its individual
!  constituent chemistry species before each SMVGEAR chemistry timestep.
!  (bdf, bmy, 4/1/03, 1/7/09)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) NTRACER (INTEGER) : Number of tracers
!  (2 ) STT     (TYPE (XPLEX) ) : Tracer concentrations [kg/box]
!  (3 ) XNUMOL  (TYPE (XPLEX) ) : Array of molecules tracer / kg tracer 
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) STT     (TYPE (XPLEX) ) : Updated tracer concentrations [molec/cm3/box]
!
!  NOTES:
!  (1 ) Now make CSAVE a local dynamic array.  Updated comments, cosmetic 
!        changes (bmy, 4/24/03)
!  (2 ) Add OpenMP parallelization commands (bmy, 8/1/03)
!  (3 ) Now dimension args XNUMOL, STT w/ NTRACER and not NNPAR (bmy, 7/20/04)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Resize CSAVE to save local memory, for SUN compiler. (bmy, 7/14/06)
!  (6 ) Now do safe division to eliminate FP errors (phs, bmy, 2/26/08)
!  (7 ) Now change error stop 30000 into a warning (phs, ccc, bmy, 1/7/09)
!  (8 ) Add support for adjoint calculation.  Save partitioning decision in 
!        PART_CASE:
!          = 1  ... partitioned NOX first 
!          = 2  ... partitioned OX first
!        (dkh, 07/22/05, dkh, 07/31/09) 
!
!******************************************************************************
!
      ! References to F90 modules 
      USE COMODE_MOD,   ONLY : CSPEC,     JLOP,       VOLUME
      USE ERROR_MOD,    ONLY : ALLOC_ERR, ERROR_STOP, SAFE_DIV
      USE TRACERID_MOD, ONLY : IDTOX,     IDTNOX,     IDTRMB
      USE TRACERID_MOD, ONLY : IDO3,      IDNO,       IDHNO2
      USE TRACERID_MOD, ONLY : CTRMB,     NMEMBER
      ! adj_group
      USE CHECKPT_MOD,  ONLY    : PART_CASE    
      USE LOGICAL_ADJ_MOD, ONLY : LADJ
      ! dkh debug
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD, LFD 
      USE LOGICAL_ADJ_MOD, ONLY : LPRINTFD
      

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "comode.h"

      ! Arguments
      INTEGER, INTENT(IN)    :: NTRACER
      TYPE (XPLEX),  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,NTRACER)
      TYPE (XPLEX),  INTENT(IN)    :: XNUMOL(NTRACER)

      ! Local variables
      INTEGER                :: I, J, L, N, JLOOP, IPL, JJ, KK
      INTEGER                :: CSAVEID(IGAS)
      INTEGER                :: CSAVEID_JJ(IGAS)
      INTEGER                :: CS, IDNUM, AS  
      TYPE (XPLEX)                 :: CONCTMP, CONCNOX, SUMM, SUMM1
      TYPE (XPLEX)                 :: CSAVE( ITLOOP, NTRACER )
      TYPE (XPLEX)                 :: QTEMP

      !=================================================================
      ! PARTITION begins here!
      !
      ! Copy values of CSPEC that need to be saved  (bdf, 3/30/99)
      !=================================================================

      ! Initialize
      IDNUM         = 0
      CSAVEID(:)    = 0
      CSAVEID_JJ(:) = 0

      ! Loop over tracers
      DO N = 1, NTRACER

         ! Skip if this is not a valid tracer
         IF ( IDTRMB(N,1) == 0 ) CYCLE

         ! Handle all other tracers except Ox 
         IF ( N /= IDTOX ) THEN
            DO KK = 1, NMEMBER(N)
               IDNUM             = IDNUM + 1
               JJ                = IDTRMB(N,KK)
               CSAVEID(JJ)       = IDNUM
               CSAVEID_JJ(IDNUM) = JJ
            ENDDO

         ! Handle Ox
         ELSE IF ( IDTOX /= 0 ) THEN
            JJ                = IDTRMB(N,1)
            IDNUM             = IDNUM + 1
            CSAVEID(JJ)       = IDNUM
            CSAVEID_JJ(IDNUM) = JJ
         ENDIF
      ENDDO

      ! Loop over tracer members and boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, JLOOP )
!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1, IDNUM
      DO L = 1, NPVERT
      DO J = 1, NLAT
      DO I = 1, NLONG

         ! 1-D SMVGEAR grid box index
         JLOOP = JLOP(I,J,L)
         IF ( JLOOP == 0 ) CYCLE

         ! Store into CSAVE
         CSAVE(JLOOP,N) = CSPEC(JLOOP,CSAVEID_JJ(N))
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! dkh debug
      !IF ( LPRINTFD ) THEN 
      IF ( LPRINTFD .and. JLOP(IFD,JFD,LFD) > 0 ) THEN
         print*, ' CSPEC in partition = ', CSPEC(JLOP(IFD,JFD,LFD),:)
         print*, ' STT   in partition = ', STT(IFD,JFD,LFD,:)
         print*, ' JLOP  in partition = ', JLOP(IFD,JFD,LFD)
      ENDIF 

      !=================================================================
      ! Split each tracer up into its components (if any)
      ! Family tracers are partitioned among members according to 
      ! initial ratios. In tracer sequence, OX must be after NOX, 
      ! otherwise, adjust the code
      !=================================================================
      DO N = 1, NTRACER

         ! Skip if it's not a valid tracer
         IF ( IDTRMB(N,1) == 0 ) CYCLE

         !### Debug
         !WRITE(6,*) 'IN PARTITION N= ', N

         ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, JLOOP, CONCTMP, SUMM, KK, JJ, SUMM1, CONCNOX )
!$OMP+PRIVATE( QTEMP )
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1, NPVERT
         DO J = 1, NLAT
         DO I = 1, NLONG

            ! 1-D SMVGEAR grid box index
            JLOOP = JLOP(I,J,L)
            IF ( JLOOP == 0 ) CYCLE

            ! Convert tracer concentration from [kg/box] to [molec/cm3/box]
            STT(I,J,L,N) = STT(I,J,L,N) / VOLUME(JLOOP) * XNUMOL(N)

            ! Store concentration of tracer N at grid box (I,J,L) in CONCTMP 
            CONCTMP = STT(I,J,L,N)

            !===========================================================
            ! First, find sum of starting concentrations
            !===========================================================

            !------------------------
            ! All tracers except Ox
            !------------------------
            IF ( N /= IDTOX ) THEN
               SUMM = 0.d0

               DO KK = 1, NMEMBER(N)
                  JJ = IDTRMB(N, KK)

                  ! Error check
                  IF ( JJ == 0 ) THEN
!$OMP CRITICAL
                     PRINT *,JJ,JLOOP,N,KK,IDTRMB(N, KK)
!$OMP END CRITICAL
                  ENDIF

                  SUMM = SUMM+CSAVE(JLOOP,CSAVEID(JJ)) * (CTRMB(N,KK)+1)
               ENDDO

            !------------------------
            ! Ox
            !------------------------
            ELSE IF ( IDTOX /= 0 ) THEN
               JJ   = IDTRMB(N,1)
               SUMM  = CSAVE(JLOOP,CSAVEID(JJ)) * (CTRMB(N,1)+1)
               SUMM1 = 0.d0


               ! SUMM  = sum of starting values for all Ox species (incl. O3)
               ! SUMM1 = sum of new values for all Ox species except O3,
               ! based on NOx partitioning
               DO KK = 2, NMEMBER(N)
                  JJ   = IDTRMB(N,KK)
                  SUMM  = SUMM+ CSAVE(JLOOP,CSAVEID(JJ))*(CTRMB(N,KK)+1)
                  SUMM1 = SUMM1+ CSPEC(JLOOP,JJ) * (CTRMB(N,KK)+1)
               ENDDO

            ENDIF

            !===========================================================
            ! Now perform the partitioning
            !===========================================================

            !----------------------------------------
            ! All tracers except Ox
            !----------------------------------------
            IF ( N /= IDTOX ) THEN

               ! Loop over # of member species in this tracer
               DO KK = 1, NMEMBER(N)

                  ! Index of member species for CSPEC
                  JJ = IDTRMB(N, KK)

                  ! QTEMP is the fraction of the given member species KK
                  ! in the tracer N.  The value QTEMP*CONCTMP is the 
                  ! concentration of the member species itself, and that
                  ! needs to be saved into CSPEC.
                  !
                  ! In the partitioning, now be sure to perform a safe
                  ! floating point division of CSAVE/SUMM.  Return the value
                  ! 1/NMEMBER(N) if the division can't be done, i.e. do a
                  ! uniform paritioning among all member species of the
                  ! given tracer. (phs, bmy, 2/26/08)
                  QTEMP = SAFE_DIV( CSAVE(JLOOP,CSAVEID(JJ)),
     &                              SUMM, xplx(1d0/NMEMBER(N)) )
                     
                  ! Store the concentration of member species KK
                  ! into the CSPEC array.  Do not allow underflow!
                  CSPEC(JLOOP,JJ) = MAX( QTEMP*CONCTMP, SMAL2 )
               ENDDO

            !----------------------------------------
            ! For Ox, take O3 = Ox - SUMM(NO2+NO3*2)
            !----------------------------------------
            ELSE IF ( IDTOX /= 0 .AND. IDTNOX /= 0 ) THEN

               ! Find Ox in CSPEC
               JJ              = IDO3
               CSPEC(JLOOP,JJ) = CONCTMP - SUMM1

               ! If Ox partitioning is OK, then skip to next tracer
               ! Old code:
               !----------------------------------------
               !IF ( CSPEC(JLOOP,JJ) > 0.0d0 ) GOTO 220
               !----------------------------------------
               ! New code: Now store PART_CASE for checkpointing. (dkh, 07/22/05)
               ! adj_group: add this to GCv8 (dkh, 07/31/09) 
               IF ( CSPEC(JLOOP,JJ) > 0.0d0 ) THEN
                  IF ( LADJ ) PART_CASE(JLOOP) = 1
                  GOTO 220
               ELSE
                  IF ( LADJ ) PART_CASE(JLOOP) = 2
               ENDIF
               !----------------------------------------

               !---------------------------------------------------------
               ! Ox partitioning failed, we are getting a negative ozone 
               ! concentration.  Instead, try partitioning Ox before NOx
               !---------------------------------------------------------

               ! Loop over member species in Ox
               DO KK = 1, NMEMBER(N)

                  ! Index of member species for CSPEC array
                  JJ = IDTRMB(N, KK)

                  ! QTEMP is the fraction of the given member species in the
                  ! Ox tracer.  The value QTEMP*CONCTMP is the concentration 
                  ! of the member species itself, and that needs to be
                  ! saved into CSPEC.
                  !
                  ! In the partitioning, now be sure to perform a safe
                  ! floating point division of CSAVE/SUMM.  Return the value
                  ! 1/NMEMBER(N) if the division can't be done, i.e. do a
                  ! uniform paritioning among all member species of the
                  ! given tracer. (phs, bmy, 2/26/08)
                  QTEMP = SAFE_DIV( CSAVE(JLOOP,CSAVEID(JJ)),
     &                              SUMM, xplx(1d0/NMEMBER(N)) )

                  ! Store the concentration of member species KK
                  ! into the CSPEC array.  Do not allow underflow!
                  CSPEC(JLOOP,JJ) = MAX( QTEMP*CONCTMP, SMAL2 )
               ENDDO

               !---------------------------------------------------------
               ! then partition NO+HNO2 
               ! (the only NOx species not contained in Ox)
               ! SUMM  = sum of starting values for NO and HNO2
               ! SUMM1 = sum of new values for all NOx species except 
               ! NO and HNO2, based on Ox partitioning
               !---------------------------------------------------------
               SUMM  = 0.d0
               SUMM1 = 0.d0
               
               ! Loop over member species of NOx
               DO KK = 1, NMEMBER(IDTNOX)
                  JJ = IDTRMB(IDTNOX, KK)

                  IF ( JJ == IDNO .OR. JJ == IDHNO2 ) THEN
                     SUMM = SUMM + CSAVE(JLOOP,CSAVEID(JJ)) *
     &                           (CTRMB(IDTNOX,KK)+1)
                  ELSE
                     SUMM1 = SUMM1+ CSPEC(JLOOP,JJ)*(CTRMB(IDTNOX,KK)+1)
                  ENDIF
               ENDDO

               ! Get NOx concentration from STT
               CONCNOX = STT(I,J,L,IDTNOX)

               ! Error test
               IF ( CONCNOX - SUMM1 < 0.d0 ) THEN
                  !------------------------------------------------------
                  ! Prior to 1/7/09
                  ! Don't stop w/ error, but just print warning msg.
                  ! Sometimes the new TPCORE can cause this error to 
                  ! trap if there CONCNOX = 0, but that can be purely 
                  ! a numerical condition and not really an error. 
                  ! (phs, ccc, bmy, 1/7/09)
                  !CALL ERROR_STOP( 'STOP 30000', 'partition.f' )
                  !------------------------------------------------------
!$OMP CRITICAL 
                  PRINT*, '### In partition.f: CONCNOX - SUMM1 < 0'
                 PRINT*, '### If CONCNOX = 0 and SUMM1 ~ 1e-99 it is OK'
                  PRINT*, '### I, J, L : ', I, J, L
                  PRINT*, '### CONCNOX : ', CONCNOX
                  PRINT*, '### SUMM1    : ', SUMM1
!$OMP END CRITICAL
               ENDIF

               ! Loop over member species in NOx
               DO KK = 1, NMEMBER(IDTNOX)

                  ! Index of member species for CSPEC
                  JJ = IDTRMB(IDTNOX,KK)

                  ! For species NO and NO2 ...
                  IF ( JJ == IDNO .OR. JJ == IDHNO2 ) THEN

                     ! QTEMP is the fraction of the given member species in 
                     ! the Ox tracer.  The value QTEMP*CONCTMP is the 
                     ! concentration of the member species itself, and that 
                     ! needs to be saved into CSPEC.
                     !
                     ! In the partitioning, now be sure to perform a safe
                     ! floating point division of CSAVE/SUMM.  Return the value
                     ! 1/NMEMBER(N) if the division can't be done, i.e. do a
                     ! uniform paritioning among all member species of the
                     ! given tracer. (phs, bmy, 2/26/08)
                     QTEMP = SAFE_DIV( CSAVE(JLOOP,CSAVEID(JJ)),
     &                               SUMM, xplx(1d0/NMEMBER(IDTNOX)))

                     ! Store the concentration of member species NO or HNO2
                     ! into the CSPEC array.  Do not allow underflow!
                     CSPEC(JLOOP,JJ) = MAX(QTEMP*(CONCNOX-SUMM1), SMAL2)
                  ENDIF
               ENDDO

               !========================================================
               ! Ox partitioning is OK
               !========================================================
 220           CONTINUE
            ENDIF
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDDO

      ! Return to calling program
      END SUBROUTINE PARTITION
