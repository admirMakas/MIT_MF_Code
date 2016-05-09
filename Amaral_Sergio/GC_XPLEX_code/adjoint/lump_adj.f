!$Id: lump_adj.f,v 1.1 2009/08/17 03:59:52 daven Exp $

      SUBROUTINE LUMP_ADJ( NTRACER, XNUMOL, STT_ADJ )
!
!******************************************************************************
!  Subroutine LUMP_ADJ takes adjoints of tracerst (STT_ADJ) and partitions them
!  into adjoints of individual chemical species (CSPEC_ADJ).  Based on 
!  ADJ_LUMP from the GCv6 adjoint (dkh, 07/31/09).
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) NTRACER (INTEGER) : Number of tracers
!  (2 ) XNUMOL  (TYPE (XPLEX) ) : Array of molecules tracer / kg tracer 
!  (3 ) STT_ADJ (TYPE (XPLEX) ) : Adjoint Tracer concentrations 
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) STT_ADJ (TYPE (XPLEX) ) : Tracer concentrations [kg/box]
!
!  Module variables included via USE as Input / Output:
!  ============================================================================
!  (1 ) CSPEC_ADJ (TYPE (XPLEX))  : Adjoint species concentrations
!
!  NOTES:
!  (1 ) Disable OMP parallel loops, which were leading to small erros in 
!       the 7th digit. (dkh, 10/08/06)  
!  (2 ) Update for GCv8 (dkh, 07/31/09) 
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,       ONLY : CSPEC,  JLOP,    VOLUME, CSPEC_ADJ
      USE TRACERID_MOD,     ONLY : IDTRMB, NMEMBER, CTRMB

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters
#     include "comode.h"     ! SMVGEAR II arrays

      ! Arguments
      INTEGER, INTENT(IN)    :: NTRACER
      TYPE (XPLEX),  INTENT(IN)    :: XNUMOL(NNPAR)
      TYPE (XPLEX),  INTENT(INOUT) :: STT_ADJ(IIPAR,JJPAR,LLPAR,NTRACER)
      ! make this an allocatable array in comode_mod
      !TYPE (XPLEX),  INTENT(INOUT) :: CSPEC_ADJ(ITLOOP,IGAS) ?

      ! Local variables
      INTEGER                :: I, J, L, N, JLOOP, KK, JJ, NN
      TYPE (XPLEX)                 :: ADCONCTMP

      !=================================================================
      ! LUMP_ADJ begins here!
      !=================================================================
      ! note: CSPEC_ADJ is initialized to zero when it is allocated.
      ! After the first call to PARTITION_ADJ it will no longer be zero 
      ! before this routine.  

!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L, N, JLOOP, ADCONCTMP, KK, JJ, NN )
!!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1, NTRACER

      ! Get index of adj species array from index of fwd species array. 
      !NN = ADJ2STT(N)
      NN = N

      ! Skip if not a valid tracer.
      ! IDTRMB for the fwd tracer (is this BCPI, etc?), NN for the adjoint tracer. 
      IF ( IDTRMB(N,1) == 0 .OR. NN == 0 ) CYCLE 
      
      DO L = 1, NPVERT
      DO J = 1, NLAT 
      DO I = 1, NLONG

         ! Initialize
         ADCONCTMP = 0.D0
 
         ! Get vector index from 3-D array indicies
         JLOOP = JLOP(I,J,L)
         IF ( JLOOP == 0 ) CYCLE

         ! Adjoint of unit conversion ( molec/cm3/box to kg/box )
         STT_ADJ(I,J,L,NN) = STT_ADJ(I,J,L,NN) * VOLUME(JLOOP) 
     &                     / XNUMOL(N)

         ADCONCTMP = ADCONCTMP + STT_ADJ(I,J,L,NN)

         ! Reset STT_ADJ to zero.  This way it won't intefere in ADJ_PARTITION
         STT_ADJ(I,J,L,NN) = 0.d0

         ! Lump adjoint values together according to families. 
         DO KK = 1, NMEMBER(N)
            JJ = IDTRMB(N,KK)
            CSPEC_ADJ(JLOOP,JJ) = CSPEC_ADJ(JLOOP,JJ) 
     &                          + ADCONCTMP * ( 1 + CTRMB(N,KK) )

         ENDDO

         ADCONCTMP = 0.D0 

      ENDDO
      ENDDO
      ENDDO
      ENDDO
!!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE LUMP_ADJ

