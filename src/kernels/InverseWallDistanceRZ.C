/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "InverseWallDistanceRZ.h"

template<>
InputParameters validParams<InverseWallDistanceRZ>()
{
  InputParameters params = validParams<Kernel>();
  params.addParam<Real>("sigma", 0.2, "smoothing parameter");
  return params;
}

InverseWallDistanceRZ::InverseWallDistanceRZ(const InputParameters & parameters) :
    Kernel(parameters),
    _sigma(getParam<Real>("sigma"))
{
}

Real
InverseWallDistanceRZ::computeQpResidual()
{
  Real r = _q_point[_qp](0);
  Real gamma = 1.0 + 2.0 * _sigma;
  return (1.0 - _sigma) * _grad_u[_qp] * _grad_u[_qp] * _test[_i][_qp] - _sigma * _u[_qp] * _grad_u[_qp] * _grad_test[_i][_qp] - gamma * std::pow(_u[_qp], 4.0) * _test[_i][_qp] +
         _sigma * ((_u[_qp] * _grad_u[_qp](0)) / r) * _test[_i][_qp]; // extra term for RZ
}

Real
InverseWallDistanceRZ::computeQpJacobian()
{
  Real r = _q_point[_qp](0);
  Real gamma = 1.0 + 2.0 * _sigma;
  return (1.0 - _sigma) * 2.0 * _grad_u[_qp] * _grad_phi[_j][_qp] * _test[_i][_qp] - _sigma * _u[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp] - _sigma * _phi[_j][_qp] * _grad_u[_qp] * _grad_test[_i][_qp] - 4.0 * gamma * std::pow(_u[_qp], 3.0) * _phi[_j][_qp] * _test[_i][_qp] +
         _sigma * ((_phi[_j][_qp] * _grad_u[_qp](0)) / r) * _test[_i][_qp] + _sigma * ((_u[_qp] * _grad_phi[_j][_qp](0)) / r) * _test[_i][_qp]; // extra terms for RZ
}
