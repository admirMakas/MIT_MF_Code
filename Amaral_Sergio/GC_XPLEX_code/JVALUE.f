C $Id: JVALUE.f,v 1.1.1.1 2009/06/09 21:51:51 daven Exp $
      SUBROUTINE JVALUE( SA )
C-----------------------------------------------------------------------
c  Calculate the actinic flux at each level for the current SZA value.
C        quit when SZA > 98.0 deg ==> tangent height = 63 km
C             or         99.                           80 km
C-----------------------------------------------------------------------
C  Add the following input variables for CTM interface (bmy, 9/13/99)
C
C  Variable  Type    Dimensn Units   Description
C  --------  ----    ------- -----   -----------
C  SA        dble    -       -       Surface Albedo
C-----------------------------------------------------------------------
c
c     AVGF   Attenuation of beam at each level for each wavelength
c     FFF    Actinic flux at each desired level
c     WAVE   Effective wavelength of each wavelength bin
c     XQO2   Absorption cross-section of O2
c     XQO3   Absorption cross-section of O3
c
C-----------------------------------------------------------------------
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "cmn_fj.h"
#     include "jv_cmn.h"

C=============== INPUT PARAMETERS ======================================
      TYPE (XPLEX), INTENT(IN) :: SA

C=============== LOCAL VARIABLES =======================================
      integer j, k
      TYPE (XPLEX)  wave, xseco3, xseco2
      TYPE (XPLEX)  AVGF(lpar),XQO3(NB),XQO2(NB)
C
      do J=1,jpnl
        do K=NW1,NW2
          FFF(K,J) = 0.d0
        enddo
      enddo
c
c---SZA check
c      write(6,1000) SZA, RFLECT, (OD(nslon,nslat,j),j=1,lpar)
      if(SZA.gt.szamax) GOTO 99
c
C---Calculate spherical weighting functions
      CALL SPHERE
c
C---Loop over all wavelength bins
      do K=NW1,NW2
        WAVE = WL(K)
        do J=1,NB
          XQO3(J) = XSECO3(K,xplx(TJ(J)))
        enddo
        do J=1,NB
          XQO2(J) = XSECO2(K,xplx(TJ(J)))
        enddo
C-----------------------------------------
        CALL OPMIE(K,WAVE,XQO2,XQO3,AVGF)
C-----------------------------------------
        do J=1,jpnl
          FFF(K,J) = FFF(K,J) + FL(K)*AVGF(J)
          !print*,'AVGF(J),FL(K),FFF(K,J)',AVGF(J),FL(K),FFF(K,J)
        enddo
      enddo
c
   99 continue
 1000 format('  SZA=',f6.1,' Reflectvty=',f6.3,' OD=',10(1pe10.3))
      return
      end
