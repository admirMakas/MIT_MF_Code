!$Id: partition_adj.f,v 1.3 2012/05/09 22:31:56 nicolas Exp $
! 
      SUBROUTINE PARTITION_ADJ( STT_ADJ, STT, NTRACER, XNUMOL ) 
!
!******************************************************************************
!  Subroutine PARTITION_ADJ is the adjoint of the fwd routine PARTITION. 
!  (dkh, 08/01/05)
!  Based on ADJ_PARTITION from the GCv6 adjoint (dkh, 07/31/09) 
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) STT     (TYPE (XPLEX) ) : Tracer concentrations [kg/box]
!  (2 ) NTRACER (INTEGER) : Number of tracers
!  (3 ) XNUMOL  (TYPE (XPLEX) ) : Array of molecules tracer / kg tracer 
!  (4 ) STT_ADJ (TYPE (XPLEX) ) : Array of adjoint concentrations
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) STT_ADJ     (TYPE (XPLEX) ) : Updated adjoint concentrations 
!
!  NOTES:
!  (1 ) See fwd version for additional notes. 
!  (2 ) Disable OMP parallel loops, which were leading to small errors
!       in the 7th digit. (dkh, 10/08/06)  
!  (3 ) Update to GCv8 (dkh, 07/31/09) 
!  (4 ) Tighten filter to 1d-10 (jkoo, dkh, boun, 05/08/12) 
!******************************************************************************
!
      ! References to F90 modules 
      USE COMODE_MOD,  ONLY : CSPEC,       JLOP,     VOLUME
      USE COMODE_MOD,  ONLY : CSPEC_PRIOR, CSPEC_ADJ
      USE ERROR_MOD,   ONLY : ALLOC_ERR,   ERROR_STOP
      USE CHECKPT_MOD, ONLY : PART_CASE
      USE TRACERID_MOD

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "comode.h"

      ! Arguments
      INTEGER, INTENT(IN)    :: NTRACER
      TYPE (XPLEX),  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,NNPAR)
      TYPE (XPLEX),  INTENT(INOUT) :: STT_ADJ(IIPAR,JJPAR,LLPAR,NTRACER)
      TYPE (XPLEX),  INTENT(IN)    :: XNUMOL(NNPAR)

      ! Local variables
      INTEGER                :: I, J, L, N, JLOOP, IPL, JJ, KK
      INTEGER                :: CSAVEID(IGAS)
      INTEGER                :: CSAVEID_JJ(IGAS)
      INTEGER                :: CS, IDNUM, AS  
      TYPE (XPLEX)                 :: CONCTMP, CONCNOX, SUMM, SUMM1
      TYPE (XPLEX)                 :: CSAVE( ITLOOP, IGAS )
      ! UPDATE: Add this so we don't overwrite CSPEC. (dkh, 10/10/08) 
      TYPE (XPLEX)                 :: CSPEC_TMP(IGAS)


      ! Adjoint varialbes 
      INTEGER :: NN
      TYPE (XPLEX)  :: ADCONCNOX
      TYPE (XPLEX)  :: ADCONCTMP
      TYPE (XPLEX)  :: ADSUMM, ADSUMM1
      TYPE (XPLEX)  :: ADCSAVE( ITLOOP, IGAS )

      
      !=================================================================
      ! PARTITION_ADJ begins here!
      !=================================================================

      ! Move this to further down below so that it happens every time
      ! through the loop. (dkh, 10/10/08) 
      !! Reset local adjoint variables to zero
      !ADSUMM     = 0.d0
      !ADSUMM1    = 0.d0
      !ADCONCNOX = 0.d0
      !ADCONCTMP = 0.d0

      ADCSAVE(:,:) = 0.d0


      ! Copy values of CSPEC that need to be saved  (bdf, 3/30/99)
      !=================================================================
      IDNUM = 0

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

      !=================================================================
      ! Split each tracer up into its components (if any)
      ! Family tracers are partitioned among members according to 
      ! initial ratios. In tracer sequence, OX must be after NOX, 
      ! otherwise, adjust the code
      !=================================================================

      ! BUG FIX: Loop from NTRACER to 1 by -1. (dkh, 10/10/08) 
      ! OLD CODE:
      !DO N = 1, NTRACER
      ! NEW CODE:
      DO N = NTRACER, 1, -1

         ! Get STT_ADJ tracer ID
         !NN = ADJ2STT(N)
         NN = N

         ! Skip if it's not a valid tracer
         IF ( IDTRMB(N,1) == 0 .OR. NN == 0 ) CYCLE

         !### Debug
         !WRITE(6,*) 'IN PARTITION N= ', N

         ! Loop over grid boxes
! UPDATE: reinstate OMP parallelization here (dkh, 10/11/08) 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, JLOOP, CONCTMP, SUMM, KK, JJ, SUMM1, CONCNOX )
!$OMP+PRIVATE( ADCONCTMP, ADCONCNOX, ADSUMM, ADSUMM1 )
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1, NPVERT
         DO J = 1, NLAT
         DO I = 1, NLONG

            ! BUG FIX: Reset adjoint variables here. (dkh, 10/10/08) 
            ! Reset local adjoint variables to zero
            ADSUMM     = 0.d0
            ADSUMM1    = 0.d0
            ADCONCNOX = 0.d0
            ADCONCTMP = 0.d0

            ! 1-D SMVGEAR grid box index
            JLOOP = JLOP(I,J,L)
            IF ( JLOOP == 0 ) CYCLE

            ! Update:  don't overwrite STT here  (dkh, 10/10/08) 
            ! OLD:
            ! Convert tracer concentration from [kg/box] to [molec/cm3/box]
            !STT(I,J,L,N) = STT(I,J,L,N) / VOLUME(JLOOP) * XNUMOL(N)
            !
            ! Store concentration in CONCTMP 
            !CONCTMP = STT(I,J,L,N)
            ! NEW:
            CONCTMP = STT(I,J,L,N) / VOLUME(JLOOP) * XNUMOL(N)

        
            ! Adjoint depends on which family was partitioned first, Ox or NOx     
            SELECT CASE ( PART_CASE(JLOOP) )

               ! Partition NOx first
              CASE ( 1 )

               !===========================================================
               ! First, find SUMM of starting concentrations
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
                        PRINT *,JJ,JLOOP,N,KK,IDTRMB(N, KK)
                     ENDIF
   
                     SUMM = SUMM 
     &                   + CSAVE(JLOOP,CSAVEID(JJ)) * (CTRMB(N,KK)+1)
                  ENDDO
               ENDIF

               ! Begin TAMC generated adjoint code.  Manual modifications
               ! are capitalized. 
               if (n .ne. idtox) then
                   do kk = 1, nmember(n)
                    ! Update: Avoid divide by SUMM**2, 
                    ! which leads to NaNs (dkh, 10/10/08) 
                    IF ( SUMM .GT. 1d-10 ) THEN
                     jj = idtrmb(n,kk)
                     adconctmp = adconctmp+cspec_adj(jloop,jj)*
     $(csave(jloop,csaveid(jj))/SUMM)
                     adcsave(jloop,csaveid(jj)) = adcsave(jloop,
     $csaveid(jj))+cspec_adj(jloop,jj)/SUMM*conctmp
                     adSUMM = adSUMM-cspec_adj(jloop,jj)*csave(jloop,
     $csaveid(jj))/(SUMM*SUMM)*conctmp
                    ENDIF
                     cspec_adj(jloop,jj) = 0.
                   end do
               else if (idtox .ne. 0 .and. idtnox .ne. 0) then
                   jj = ido3
                   adconctmp = adconctmp+cspec_adj(jloop,jj)
                   adSUMM1 = adSUMM1-cspec_adj(jloop,jj)
                   cspec_adj(jloop,jj) = 0.
               else
                  print*, ' big error here ' 
               endif
               if (n .ne. idtox) then
                   do kk = 1, nmember(n)
                     jj = idtrmb(n,kk)
                     adcsave(jloop,csaveid(jj)) = adcsave(jloop,
     $csaveid(jj))+adSUMM*(1+ctrmb(n,kk))
                   end do
                   adSUMM = 0.
               else if (idtox .ne. 0) then
                   do kk = 2, nmember(n)
                     jj = idtrmb(n,kk)
                     cspec_adj(jloop,jj) = cspec_adj(jloop,jj)
     &+adSUMM1*(1+ctrmb(n,kk))
                   end do
                   adSUMM1 = 0.
               endif

              ! Partition Ox first
              CASE ( 2 )

               !===========================================================
               ! First, find SUMM of starting concentrations
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
                        PRINT *,JJ,JLOOP,N,KK,IDTRMB(N, KK)
                     ENDIF
   
                     SUMM = SUMM 
     &                   + CSAVE(JLOOP,CSAVEID(JJ)) * (CTRMB(N,KK)+1)
                  ENDDO

               !------------------------
               ! Ox
               !------------------------
               ELSE IF ( IDTOX /= 0 ) THEN
                  JJ   = IDTRMB(N,1)
                  SUMM  = CSAVE(JLOOP,CSAVEID(JJ)) * (CTRMB(N,1)+1)
! Case 1 stuff. dkh
!               SUMM1 = 0.d0


                  ! SUMM  = SUMM of starting values for all Ox species (incl. O3)
                  ! SUMM1 = SUMM of new values for all Ox species except O3,
                  ! based on NOx partitioning
                  DO KK = 2, NMEMBER(N)
                     JJ   = IDTRMB(N,KK)
                     SUMM  = SUMM 
     &                    + CSAVE(JLOOP,CSAVEID(JJ))*(CTRMB(N,KK)+1)
! Case 1 stuff. dkh
!                  SUMM1 = SUMM1+ CSPEC(JLOOP,JJ) * (CTRMB(N,KK)+1)
                  ENDDO

               ENDIF
            
               ! Begin TAMC generated adjoint of partioning for Case 2
               if (n .ne. idtox) then
                     do kk = 1, nmember(n)
                      IF ( SUMM .gt. 1d-10 ) THEN
                       jj = idtrmb(n,kk)
                       adconctmp = adconctmp+cspec_adj(jloop,jj)*
     $(csave(jloop,csaveid(jj))/SUMM)
                       adcsave(jloop,csaveid(jj)) = adcsave(jloop,
     $csaveid(jj))+cspec_adj(jloop,jj)/SUMM*conctmp
                       adSUMM = adSUMM-cspec_adj(jloop,jj)*csave(jloop,
     $csaveid(jj))/(SUMM*SUMM)*conctmp
                      ENDIF
                       cspec_adj(jloop,jj) = 0.
                  end do
               else if (idtox .ne. 0 .and. idtnox .ne. 0) then
                     do kk = 1, nmember(n)
                       jj = idtrmb(n,kk)
                       ! Update: don't overwrite CSPEC (dkh, 10/10/08) 
                       !cspec(jloop,jj) = csave(jloop,csaveid(jj))/SUMM*
                       CSPEC_TMP(jj) = csave(jloop,csaveid(jj))/SUMM*
     $conctmp
                     end do
                     SUMM = 0.d0
                     SUMM1 = 0.d0
                     do kk = 1, nmember(idtnox)
                       jj = idtrmb(idtnox,kk)
                       if (jj .eq. idno .or. jj .eq. idhno2) then
      SUMM = SUMM+csave(jloop,csaveid(jj))*(ctrmb(idtnox,kk)+1)
                       else
      ! Update: use CSPEC_TMP (dkh, 10/10/08) 
      !SUMM1 = SUMM1+cspec(jloop,jj)*(ctrmb(idtnox,kk)+1)
      SUMM1 = SUMM1+CSPEC_TMP(jj)*(ctrmb(idtnox,kk)+1)
                       endif
                     end do
                     !---------------------------------- 
                     ! BUG FIX: need to convert units
                     ! of concnox here (jkoo, dkh, 09/30/10) 
                     ! old code:
                     !concnox = stt(i,j,l,idtnox)
                     ! new code: 
                     concnox = stt(i,j,l,idtnox)
     &                 / VOLUME(JLOOP) * XNUMOL(IDTNOX)
                     !---------------------------------- 
                     do kk = 1, nmember(idtnox)
                       jj = idtrmb(idtnox,kk)
                       if (jj .eq. idno .or. jj .eq. idhno2) then
                        !IF ( SUMM .gt. 0d0 ) THEN
                        IF ( SUMM .gt. 1d-10 ) THEN
      adconcnox = adconcnox+cspec_adj(jloop,jj)
     &*(csave(jloop,csaveid(jj))/SUMM)
      adcsave(jloop,csaveid(jj)) = adcsave(jloop,csaveid(jj))+
     $cspec_adj(jloop,jj)/SUMM*(concnox-SUMM1)
      adSUMM = adSUMM-cspec_adj(jloop,jj)*csave(jloop,csaveid(jj))
     $/(SUMM*SUMM)
     $*(concnox-SUMM1)
      adSUMM1=adSUMM1-cspec_adj(jloop,jj)*
     & (csave(jloop,csaveid(jj))/SUMM)
                        ENDIF
      cspec_adj(jloop,jj) = 0.
                       endif
                     end do
                     STT_ADJ(i,j,l,idtnox) = 
     &                STT_ADJ(i,j,l,idtnox)+adconcnox
                     adconcnox = 0.
                     do kk = 1, nmember(idtnox)
                       jj = idtrmb(idtnox,kk)
                       if (jj .eq. idno .or. jj .eq. idhno2) then
      adcsave(jloop,csaveid(jj)) = adcsave(jloop,csaveid(jj))+adSUMM*(1+
     $ctrmb(idtnox,kk))
                       else
      cspec_adj(jloop,jj) = cspec_adj(jloop,jj)
     &                    +adSUMM1*(1+ctrmb(idtnox,kk))
                       endif
                     end do
                     adSUMM1 = 0.
                     adSUMM = 0.
                     !SUMM = SUMMk
                     if (n .ne. idtox) then
                       SUMM = 0.d0
                       do kk = 1, nmember(n)
      jj = idtrmb(n,kk)
      SUMM = SUMM+csave(jloop,csaveid(jj))*(ctrmb(n,kk)+1)
                       end do
                     else if (idtox .ne. 0) then
                       jj = idtrmb(n,1)
                       SUMM = csave(jloop,csaveid(jj))*(ctrmb(n,1)+1)
                       do kk = 2, nmember(n)
      jj = idtrmb(n,kk)
      SUMM = SUMM+csave(jloop,csaveid(jj))*(ctrmb(n,kk)+1)
                       end do
                     endif
                     do kk = 1, nmember(n)
                      !IF ( SUMM .gt. 0d0 ) THEN
                      IF ( SUMM .gt. 1d-10 ) THEN
                       jj = idtrmb(n,kk)
                       adconctmp = adconctmp+cspec_adj(jloop,jj)*
     $(csave(jloop,csaveid(jj))/SUMM)
                       adcsave(jloop,csaveid(jj)) = adcsave(jloop,
     $csaveid(jj))+cspec_adj(jloop,jj)/SUMM*conctmp
                       adSUMM = adSUMM-cspec_adj(jloop,jj)*csave(jloop,
     $csaveid(jj))/(SUMM*SUMM)*conctmp
                      ENDIF
                       cspec_adj(jloop,jj) = 0.
                     end do
               endif
               if (n .ne. idtox) then
                     do kk = 1, nmember(n)
                       jj = idtrmb(n,kk)
                       adcsave(jloop,csaveid(jj)) = adcsave(jloop,
     $csaveid(jj))+adSUMM*(1+ctrmb(n,kk))
                     end do
                     adSUMM = 0.
               else if (idtox .ne. 0) then
                     do kk = 2, nmember(n)
                       jj = idtrmb(n,kk)
                       adcsave(jloop,csaveid(jj)) = adcsave(jloop,
     $csaveid(jj))+adSUMM*(1+ctrmb(n,kk))
                     end do
                     jj = idtrmb(n,1)
                     adcsave(jloop,csaveid(jj)) = 
     $adcsave(jloop,csaveid(jj))+adSUMM*(1+ctrmb(n,1))
                     adSUMM = 0.
               endif

              CASE DEFAULT
               WRITE(6,*) I, J, L, JLOOP
               CALL ERROR_STOP( 'bad PART_CASE', 'PARTITION_ADJ' )

            END SELECT
           

            STT_ADJ(i,j,l,NN) = STT_ADJ(i,j,l,NN)+adconctmp
            adconctmp = 0.
            STT_ADJ(i,j,l,NN) = STT_ADJ(i,j,l,NN)
     &                        /volume(jloop)*xnumol(n)

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDDO

      ! BUG FIX: add this part to pass sensitivities to ADJ_CSPEC 
      ! which then get fed back to STT_ADJ in subroutine LUMP_ADJ
      ! (dkh, 10/11/08) 
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
         
         ! fwd code:
         !CSAVE(JLOOP,N) = CSPEC(JLOOP,CSAVEID_JJ(N))
         ! adj code:
         CSPEC_ADJ(JLOOP,CSAVEID_JJ(N)) = CSPEC_ADJ(JLOOP,CSAVEID_JJ(N))
     &                                  + ADCSAVE(JLOOP,N) 

      ENDDO 
      ENDDO 
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program 
      END SUBROUTINE PARTITION_ADJ

