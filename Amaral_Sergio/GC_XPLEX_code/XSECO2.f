C $Id: XSECO2.f,v 1.1.1.1 2009/06/09 21:51:52 daven Exp $
      FUNCTION XSECO2(K,TTT)
C-----------------------------------------------------------------------
c  Cross-sections for O2 interpolated across 3 temps; No S_R Bands yet!
C-----------------------------------------------------------------------
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
      
#     include "cmn_fj.h"
#     include "jv_cmn.h"

      integer k
      TYPE (XPLEX) ttt, flint, xseco2
      XSECO2 =
     F  FLINT(TTT,TQQ(1,1),TQQ(2,1),TQQ(3,1),QO2(K,1),QO2(K,2),QO2(K,3))
      return
      end
