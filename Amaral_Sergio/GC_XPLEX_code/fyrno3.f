! $Id: fyrno3.f,v 1.1.1.1 2009/06/09 21:51:52 daven Exp $
      TYPE (XPLEX) FUNCTION FYRNO3( XCARBN, ZDNUM, TT )
!
!******************************************************************************
!  Function FYRNO3 returns organic nitrate yields YN = RKA/(RKA+RKB)
!  from RO2+NO reactions as a function of the number N of carbon atoms.
!  (lwh, jyl, gmg, djj, bmy, 1/1/89, 6/26/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) XCARBN (TYPE (XPLEX)) : Number of C atoms in RO2
!  (2 ) ZDNUM  (TYPE (XPLEX)) : Air density   [molec/cm3 ]
!  (3 ) TT     (TYPE (XPLEX)) : Temperature   [K         ]
!    
!  NOTES: 
!  (1 ) Original code from Larry Horowitz, Jinyou Liang, Gerry Gardner,
!        and Daniel Jacob circa 1989/1990.
!  (2 ) Updated following Atkinson 1990.
!  (3 ) Change yield from Isoprene Nitrate (ISN2) from 0.44% to 12%,
!        according to Sprengnether et al., 2002. (amf, bmy, 1/7/02)
!  (4 ) Eliminate obsolete code from 1/02 (bmy, 2/27/02)
!  (5 ) Updated comment description of XCARBN (bmy, 6/26/03)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      ! Arguments
      TYPE (XPLEX), INTENT(IN) :: XCARBN, ZDNUM, TT

      ! Local variables
      TYPE (XPLEX)             :: YYYN, XXYN,  AAA,  RARB, ZZYN
      TYPE (XPLEX)             :: XF,   ALPHA, Y300, BETA, XMINF, XM0

      ! Initialize static variables
      DATA Y300,ALPHA,BETA,XM0,XMINF,XF/
     & xplex(.826d0,0d0),xplex(1.94D-22,0d0),
     & xplex(.97d0,0d0),xplex(0.d0,0d0),xplex(8.1d0,0d0),
     & xplex(.411d0,0d0)/
      
      !=================================================================
      ! FYRNO3 begins here!
      !=================================================================
      XXYN   = ALPHA*EXP(BETA*XCARBN)*ZDNUM*((300.d0/TT)**XM0)
      YYYN   = Y300*((300.d0/TT)**XMINF)
      AAA    = LOG10(XXYN/YYYN)
      ZZYN   = 1.d0/(1.d0+ AAA*AAA )
      RARB   = (XXYN/(1.d0+ (XXYN/YYYN)))*(XF**ZZYN)
      FYRNO3 = RARB/(1.d0 + RARB)
     
      ! Return to calling program
      END FUNCTION FYRNO3
