/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "InverseWallDistanceGrad.h"

template<>
InputParameters validParams<InverseWallDistanceGrad>()
{
  InputParameters params = validParams<Kernel>();
  return params;
}

InverseWallDistanceGrad::InverseWallDistanceGrad(const InputParameters & parameters) :
    Kernel(parameters)
{
}

Real
InverseWallDistanceGrad::computeQpResidual()
{
  return _grad_u[_qp] * _grad_u[_qp] * _test[_i][_qp] - std::pow(_u[_qp], 4.0) * _test[_i][_qp];
}

Real
InverseWallDistanceGrad::computeQpJacobian()
{
  return 2.0 * _grad_u[_qp] * _grad_phi[_j][_qp] * _test[_i][_qp] - 4.0 * std::pow(_u[_qp], 3.0) * _phi[_j][_qp] * _test[_i][_qp];
}
