C $Id: GAUSSP.f,v 1.1.1.1 2009/06/09 21:51:52 daven Exp $
      SUBROUTINE GAUSSP (N,XPT,XWT)
C-----------------------------------------------------------------------
C  Loads in pre-set Gauss points for 4 angles from 0 to +1 in cos(theta)=mu
C-----------------------------------------------------------------------
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
      INTEGER N,I
      TYPE (XPLEX)  XPT(N),XWT(N)
      TYPE (XPLEX) GPT4(4),GWT4(4)
      DATA GPT4%r/.06943184420297D0,.33000947820757D0,.66999052179243D0,
     G          .93056815579703D0/
      GPT4%i = 0d0
      DATA GWT4%r/.17392742256873D0,.32607257743127D0,.32607257743127D0,
     W          .17392742256873D0/
      GWT4%i = 0d0
      N = 4
      DO I=1,N
        XPT(I) = GPT4(I)
        XWT(I) = GWT4(I)
      ENDDO
      RETURN
      END
