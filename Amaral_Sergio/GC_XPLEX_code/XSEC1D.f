C $Id: XSEC1D.f,v 1.1.1.1 2009/06/09 21:51:54 daven Exp $
      FUNCTION XSEC1D(K,TTT)
C-----------------------------------------------------------------------
c  Quantum yields for O3 --> O2 + O(1D) interpolated across 3 temps
C-----------------------------------------------------------------------
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "cmn_fj.h"
#     include "jv_cmn.h"

      integer k
      TYPE (XPLEX) ttt, flint, xsec1d
      XSEC1D =
     F  FLINT(TTT,TQQ(1,3),TQQ(2,3),TQQ(3,3),Q1D(K,1),Q1D(K,2),Q1D(K,3))
      return
      end
