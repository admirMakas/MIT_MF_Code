! $Id: lump.f,v 1.1.1.1 2009/06/09 21:51:52 daven Exp $
      SUBROUTINE LUMP( NTRACER, XNUMOL, STT )
!
!******************************************************************************
!  Subroutine LUMP takes individual chemistry species and "lumps" them back 
!  into tracers after each SMVGEAR chemistry timestep. (bmy, 4/1/03, 7/20/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) NTRACER (INTEGER) : Number of tracers
!  (2 ) XNUMOL  (TYPE (XPLEX) ) : Array of molecules tracer / kg tracer 
!  (3 ) STT     (TYPE (XPLEX) ) : Tracer concentrations [molec/cm3/box]
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) STT     (TYPE (XPLEX) ) : Tracer concentrations [kg/box]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes (bmy, 4/1/03)
!  (2 ) Added OpenMP parallelization commands (bmy, 8/1/03)
!  (3 ) Now dimension args XNUMOL, STT w/ NTRACER and not NNPAR (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,   ONLY : CSPEC,  JLOP,    VOLUME
      USE TRACERID_MOD, ONLY : IDTRMB, NMEMBER, CTRMB

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters
#     include "comode.h"     ! SMVGEAR II arrays

      ! Arguments
      INTEGER, INTENT(IN)    :: NTRACER
      TYPE (XPLEX),  INTENT(IN)    :: XNUMOL(NTRACER)
      TYPE (XPLEX),  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,NTRACER)   

      ! Local variables
      INTEGER                :: I, J, L, N, JLOOP, KK, JJ
      TYPE (XPLEX)                 :: CONCTMP  
      !DO I=1,SIZE(STT,1)
      !DO J=1,SIZE(STT,2)
      !DO L=1,SIZE(STT,3)
      !DO N=1,SIZE(STT,4)
      !IF (dimag(STT(I,J,L,N))>1d-80 ) then
      !     print*,'STT NOX BEG LUMP',STT(I,J,L,N),n,i,j,l
      !ENDIF
      !ENDDO
      !ENDDO
      !ENDDO
      !ENDDO
      !=================================================================
      ! LUMP begins here!
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, JLOOP, CONCTMP, KK, JJ )
!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1, NTRACER
         
         ! Skip if not a valid tracer
         IF ( IDTRMB(N,1) == 0 ) CYCLE
       
         ! Loop over grid boxes
         DO L = 1, NPVERT
         DO J = 1, NLAT
         DO I = 1, NLONG

            ! 1-D SMVGEAR grid box index 
            JLOOP = JLOP(I,J,L)
            IF ( JLOOP == 0 ) CYCLE
            
            ! Compute tracer concentration [molec/cm3/box] by
            ! looping over all species belonging to this tracer
            CONCTMP = 0.d0
            DO KK = 1, NMEMBER(N)
               JJ = IDTRMB(N, KK)
               CONCTMP = CONCTMP + ( 1d0+CTRMB(N,KK) ) * CSPEC(JLOOP,JJ)
               !IF (dimag(CONCTMP)>1d-80 ) then
               !print*,'IMAG CONCTMP',dimag(CONCTMP)
               !print*,'IMAG CTRMB(N,KK)',dimag(CTRMB(N,KK))
               !print*,'IMAG CSPEC',dimag(CSPEC(JLOOP,J))
               !ENDIF
            ENDDO
            !IF (dimag(STT(I,J,L,N))>1d-80 ) then
            !   print*,'STT bef CONCTMP  LUMP',STT(I,J,L,N),n,i,j,l
            !ENDIF
            ! Save tracer concentrations back to STT
            STT(I,J,L,N) = CONCTMP
            !IF (dimag(STT(I,J,L,N))>1d-80 ) then
            !   print*,'STT aft CONCTMP  LUMP',STT(I,J,L,N),n,i,j,l
            !ENDIF
            ! Change STT from [molec/cm3/box] back to [kg/box]
            STT(I,J,L,N) = STT(I,J,L,N) * VOLUME(JLOOP) / XNUMOL(N)
            !IF (dimag(STT(I,J,L,N))>1d-80 ) then
            !   print*,'STT aft change units  LUMP',STT(I,J,L,N),n,i,j,l
            !ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
      !DO I=1,SIZE(STT,1)
      !DO J=1,SIZE(STT,2)
      !DO L=1,SIZE(STT,3)
      !DO N=1,SIZE(STT,4)
      !IF (dimag(STT(I,J,L,N))>1d-80 ) then
      !     print*,'STT END LUMP',STT(I,J,L,N),n,i,j,l
      !ENDIF
      !ENDDO
      !ENDDO
      !ENDDO
      !ENDDO
      ! Return to calling program
      END SUBROUTINE LUMP

