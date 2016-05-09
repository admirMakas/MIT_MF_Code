MODULE gckpp_adj_Model

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Completely defines the model gckpp_adj
!    by using all the associated modules
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE gckpp_adj_Precision
  USE gckpp_adj_Parameters
  USE gckpp_adj_Global
  USE gckpp_adj_Function
  USE gckpp_adj_Integrator
  USE gckpp_adj_Rates
  USE gckpp_adj_Jacobian
  USE gckpp_adj_Hessian
  USE gckpp_adj_Stoichiom
  USE gckpp_adj_LinearAlgebra
  USE gckpp_adj_Monitor
  USE gckpp_adj_Util

END MODULE gckpp_adj_Model

