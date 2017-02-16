/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "INSMomentumLaplaceFormWALE.h"

template<>
InputParameters validParams<INSMomentumLaplaceFormWALE>()
{
  InputParameters params = validParams<INSMomentumBaseWALE>();
  return params;
}



INSMomentumLaplaceFormWALE::INSMomentumLaplaceFormWALE(const InputParameters & parameters) :
  INSMomentumBaseWALE(parameters)
{
}



Real INSMomentumLaplaceFormWALE::computeQpResidualViscousPart()
{
  // Simplified version: mu * Laplacian(u_component)
  Real _mu = INSMomentumBaseWALE::computeQpDynamicViscosity();
  return _mu * (_grad_u[_qp] * _grad_test[_i][_qp]);
}



Real INSMomentumLaplaceFormWALE::computeQpJacobianViscousPart()
{
  // Viscous part, Laplacian version
  Real _mu = INSMomentumBaseWALE::computeQpDynamicViscosity();
  return _mu * (_grad_phi[_j][_qp] * _grad_test[_i][_qp]);
}



Real INSMomentumLaplaceFormWALE::computeQpOffDiagJacobianViscousPart(unsigned /*jvar*/)
{
  return 0.;
}
