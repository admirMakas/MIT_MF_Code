!$Id: gckpp_adj_Precision.f90,v 1.2 2009/06/12 01:44:48 daven Exp $
MODULE gckpp_adj_Precision

!
! Definition of different levels of accuracy
! for REAL variables using KIND parameterization
!
! KPP SP - Single precision kind
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,30)
! KPP DP - TYPE (XPLEX) kind
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)
! KPP QP - Quadruple precision kind
  INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(18,400)

END MODULE gckpp_adj_Precision


