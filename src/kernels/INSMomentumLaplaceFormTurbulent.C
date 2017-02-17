/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "INSMomentumLaplaceFormTurbulent.h"

template<>
InputParameters validParams<INSMomentumLaplaceFormTurbulent>()
{
  InputParameters params = validParams<INSMomentumBaseTurbulent>();
  return params;
}



INSMomentumLaplaceFormTurbulent::INSMomentumLaplaceFormTurbulent(const InputParameters & parameters) :
  INSMomentumBaseTurbulent(parameters)
{
}



Real INSMomentumLaplaceFormTurbulent::computeQpResidualViscousPart()
{
  // Simplified version: mu * Laplacian(u_component)
  return _mu[_qp] * (_grad_u[_qp] * _grad_test[_i][_qp]);
}



Real INSMomentumLaplaceFormTurbulent::computeQpJacobianViscousPart()
{
  // Viscous part, Laplacian version
  return _mu[_qp] * (_grad_phi[_j][_qp] * _grad_test[_i][_qp]);
}



Real INSMomentumLaplaceFormTurbulent::computeQpOffDiagJacobianViscousPart(unsigned /*jvar*/)
{
  return 0.;
}
