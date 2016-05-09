C $Id: SPHERE.f,v 1.1.1.1 2009/06/09 21:51:50 daven Exp $
      SUBROUTINE SPHERE
C-----------------------------------------------------------------------
c  Calculation of spherical geometry; derive tangent heights, slant path
c  lengths and air mass factor for each layer. Not called when
c  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
c  beam (where tangent height is below altitude J-value desired at).
C-----------------------------------------------------------------------
c
c     GMU     MU, cos(solar zenith angle)
c     RZ      Distance from centre of Earth to each point (cm)
c     RQ      Square of radius ratios
c     TANHT   Tangent height for the current SZA
c     XL      Slant path between points
c     AMF     Air mass factor for slab between level and level above
c
C-----------------------------------------------------------------------
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "cmn_fj.h"
#     include "jv_cmn.h"

      integer i, j, k, ii
      TYPE (XPLEX) airmas, gmu, xmu1, xmu2, xl, diff
      TYPE (XPLEX) Ux,H,RZ(NB),RQ(NB),ZBYR
c
c  Inlined air mass factor function for top of atmosphere
!      AIRMAS(Ux,H) = (1.0d0+H)/SQRT(Ux*Ux+2.0d0*H*(1.0d0-
!     $         0.6817d0*EXP(-57.3d0*ABS(Ux)/SQRT(1.0d0+5500.d0*H))/
!     $                                             (1.0d0+0.625d0*H)))
c
      GMU = U0
      RZ(1)=RAD+Z(1)
      ZBYR = ZZHT/RAD
      DO 2 II=2,NB
        RZ(II) = RAD + Z(II)
        RQ(II-1) = (RZ(II-1)/RZ(II))**2
    2 CONTINUE
      IF (GMU.LT.0.0D0) THEN
        TANHT = RZ(nlbatm)/SQRT(1.0D0-GMU**2)
      ELSE
        TANHT = RZ(nlbatm)
      ENDIF
c
c  Go up from the surface calculating the slant paths between each level
c  and the level above, and deriving the appropriate Air Mass Factor
      DO 16 J=1,NB
        DO K=1,NB
          AMF(K,J)=0.D0
        ENDDO
c
c  Air Mass Factors all zero if below the tangent height
        IF (RZ(J).LT.TANHT) GOTO 16
c  Ascend from layer J calculating AMFs
        XMU1=ABS(GMU)
        DO 12 I=J,lpar
          XMU2=SQRT(1.0D0-RQ(I)*(1.0D0-XMU1**2))
          XL=RZ(I+1)*XMU2-RZ(I)*XMU1
          AMF(I,J)=XL/(RZ(I+1)-RZ(I))
          XMU1=XMU2
   12   CONTINUE
c  Use function and scale height to provide AMF above top of model
        AMF(NB,J)=(1.0d0+ZBYR)/SQRT(XMU1*XMU1+2.0d0*ZBYR*(1.0d0-
     $         0.6817d0*EXP(-57.3d0*ABS(XMU1)/SQRT(1.0d0+5500.d0*ZBYR))/
     $                               (1.0d0+0.625d0*ZBYR)))!AIRMAS(XMU1,ZBYR)
c
c  Twilight case - Emergent Beam
        IF (GMU.GE.0.0D0) GOTO 16
        XMU1=ABS(GMU)
c  Descend from layer J
        DO 14 II=J-1,1,-1
          DIFF=RZ(II+1)*SQRT(1.0D0-XMU1**2)-RZ(II)
          if(II.eq.1) DIFF=max(DIFF,0.d0)   ! filter
c  Tangent height below current level - beam passes through twice
          IF (DIFF.LT.0.0D0) THEN
            XMU2=SQRT(1.0D0-(1.0D0-XMU1**2)/RQ(II))
            XL=ABS(RZ(II+1)*XMU1-RZ(II)*XMU2)
            AMF(II,J)=2.d0*XL/(RZ(II+1)-RZ(II))
            XMU1=XMU2
c  Lowest level intersected by emergent beam
          ELSE
            XL=RZ(II+1)*XMU1*2.0D0
c            WTING=DIFF/(RZ(II+1)-RZ(II))
c            AMF(II,J)=(1.0D0-WTING)*2.D0**XL/(RZ(II+1)-RZ(II))
            AMF(II,J)=XL/(RZ(II+1)-RZ(II))
            GOTO 16
          ENDIF
   14   CONTINUE
c
   16 CONTINUE
      RETURN
      END
