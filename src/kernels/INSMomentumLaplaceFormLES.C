/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "INSMomentumLaplaceFormLES.h"

template<>
InputParameters validParams<INSMomentumLaplaceFormLES>()
{
  InputParameters params = validParams<INSMomentumBaseLES>();
  return params;
}



INSMomentumLaplaceFormLES::INSMomentumLaplaceFormLES(const InputParameters & parameters) :
  INSMomentumBaseLES(parameters)
{
}



Real INSMomentumLaplaceFormLES::computeQpResidualViscousPart()
{
  // Simplified version: mu * Laplacian(u_component)
  return _mu[_qp] * (_grad_u[_qp] * _grad_test[_i][_qp]);
}



Real INSMomentumLaplaceFormLES::computeQpJacobianViscousPart()
{
  // Viscous part, Laplacian version
  return _mu[_qp] * (_grad_phi[_j][_qp] * _grad_test[_i][_qp]);
}



Real INSMomentumLaplaceFormLES::computeQpOffDiagJacobianViscousPart(unsigned /*jvar*/)
{
  return 0.;
}
