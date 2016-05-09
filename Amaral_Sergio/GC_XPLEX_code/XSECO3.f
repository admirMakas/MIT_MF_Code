C $Id: XSECO3.f,v 1.1.1.1 2009/06/09 21:51:52 daven Exp $
      FUNCTION XSECO3(K,TTT)
C-----------------------------------------------------------------------
c  Cross-sections for O3 for all processes interpolated across 3 temps
C-----------------------------------------------------------------------
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "cmn_fj.h"
#     include "jv_cmn.h"

      integer k
      TYPE (XPLEX) ttt, flint, xseco3
      XSECO3  = 
     F  FLINT(TTT,TQQ(1,2),TQQ(2,2),TQQ(3,2),QO3(K,1),QO3(K,2),QO3(K,3))
      return
      end
