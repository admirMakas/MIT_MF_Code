C $Id: FLINT.f,v 1.1.1.1 2009/06/09 21:51:54 daven Exp $
      TYPE (XPLEX) FUNCTION FLINT (TINT,T1,T2,T3,F1,F2,F3)
      USE MYTYPE
      USE COMPLEXIFY
C-----------------------------------------------------------------------
c  Three-point linear interpolation function
C-----------------------------------------------------------------------
      TYPE (XPLEX) TINT,T1,T2,T3,F1,F2,F3
      IF (TINT .LE. T2)  THEN
        IF (TINT .LE. T1)  THEN
          FLINT  = F1
        ELSE
          FLINT = F1 + (F2 - F1)*(TINT -T1)/(T2 -T1)
        ENDIF
      ELSE
        IF (TINT .GE. T3)  THEN
          FLINT  = F3
        ELSE
          FLINT = F2 + (F3 - F2)*(TINT -T2)/(T3 -T2)
        ENDIF
      ENDIF
      return
      end
