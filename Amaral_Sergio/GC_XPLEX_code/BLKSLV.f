C $Id: BLKSLV.f,v 1.1.1.1 2009/06/09 21:51:52 daven Exp $
      SUBROUTINE BLKSLV
C-----------------------------------------------------------------------
C  Solves the block tri-diagonal system:
C               A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
C-----------------------------------------------------------------------
      USE MYTYPE
      USE ERROR_MOD
      USE COMPLEXIFY
      IMPLICIT NONE
#     include "jv_mie.h"
      integer i, j, k, id
      TYPE (XPLEX)  SUMM
C-----------UPPER BOUNDARY ID=1
      CALL GEN(1)
      CALL MATIN4 (B)
      DO I=1,N
         RR(I,1) = 0.0d0
        DO J=1,N
          SUMM = 0.0d0
         DO K=1,N
          SUMM = SUMM - B(I,K)*CC(K,J)
!          if (isnan(B(I,K)).or.isnan(CC(K,J))) then
!             print*,'B(I,K),C(K,J)',B(I,K),CC(K,J)
!          endif 
         ENDDO
         DD(I,J,1) = SUMM
         RR(I,1) = RR(I,1) + B(I,J)*H(J)
!          if (isnan(RR(I,1))) then
!             print*,'B(I,J),H(J) blkslv 30',B(I,J),H(J)
!             CALL GEOS_CHEM_STOP
!          endif
        ENDDO
      ENDDO
C----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
      DO ID=2,ND-1
        CALL GEN(ID)
        DO I=1,N
          DO J=1,N
          B(I,J) = B(I,J) + A(I)*DD(I,J,ID-1)
          ENDDO
          H(I) = H(I) - A(I)*RR(I,ID-1)
!          if (isnan(H(I))) then
!             print*,'A(I),RR(I,ID-1)  blkslv 49',A(I),RR(I,ID-1),H(I)
!             CALL GEOS_CHEM_STOP
!          endif
        ENDDO
        CALL MATIN4 (B)
        DO I=1,N
          RR(I,ID) = 0.0d0
          DO J=1,N
          RR(I,ID) = RR(I,ID) + B(I,J)*H(J)
!          if (isnan(RR(I,ID))) then
!             print*,'B(I,J),H(J) blkslv 49',B(I,J),H(J)
!             CALL GEOS_CHEM_STOP
!          endif
          DD(I,J,ID) = - B(I,J)*C1(J)
          ENDDO
        ENDDO
      ENDDO
C---------FINAL DEPTH POINT: ND
      CALL GEN(ND)
      DO I=1,N
        DO J=1,N
          SUMM = 0.0d0
          DO K=1,N
          SUMM = SUMM + AA(I,K)*DD(K,J,ND-1)
          ENDDO
        B(I,J) = B(I,J) + SUMM
        H(I) = H(I) - AA(I,J)*RR(J,ND-1)
        ENDDO
      ENDDO
      CALL MATIN4 (B)
      DO I=1,N
        RR(I,ND) = 0.0d0
        DO J=1,N
        RR(I,ND) = RR(I,ND) + B(I,J)*H(J)
!        if (isnan(RR(I,ND))) then
!             print*,'B(I,J),H(J) blkslv 73',B(I,J),H(J)
!             CALL GEOS_CHEM_STOP
!          endif
        ENDDO
      ENDDO
C-----------BACK SOLUTION
      DO ID=ND-1,1,-1
       DO I=1,N
        DO J=1,N
         RR(I,ID) = RR(I,ID) + DD(I,J,ID)*RR(J,ID+1)
!         if (isnan(RR(I,1))) then
!             print*,'B(I,J),H(J) blkslv 83',DD(I,J,ID),H(J)
!             CALL GEOS_CHEM_STOP
!          endif
        ENDDO
       ENDDO
      ENDDO
C----------MEAN J & H
      DO ID=1,ND,2
        FJ(ID) = 0.0d0
       DO I=1,N
        FJ(ID) = FJ(ID) + RR(I,ID)*WT(I)
       ENDDO
      ENDDO
      DO ID=2,ND,2
        FJ(ID) = 0.0d0
       DO I=1,N
        FJ(ID) = FJ(ID) + RR(I,ID)*WT(I)*EMU(I)
!        if (isnan(FJ(ID))) then
!            print*,'RR,WT,EMU',RR(I,ID),WT(I),EMU(I)
!        endif
       ENDDO
      ENDDO
C Output fluxes for testing purposes
c      CALL CH_FLUX
c
      RETURN
      END
