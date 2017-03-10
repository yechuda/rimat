/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "NeeKovasznayRZ.h"

template<>
InputParameters validParams<NeeKovasznayRZ>()
{
  InputParameters params = validParams<Kernel>();

  return params;
}



NeeKovasznayRZ::NeeKovasznayRZ(const InputParameters & parameters) :
  Kernel(parameters)

{
}



Real NeeKovasznayRZ::computeQpResidual()
{
  Real r = _q_point[_qp](0);
  return _u[_qp] * _grad_u[_qp](0) / r * _test[_i][_qp];
}



Real NeeKovasznayRZ::computeQpJacobian()
{
  Real r = _q_point[_qp](0);
  return _phi[_j][_qp] * _grad_phi[_j][_qp](0) / r * _test[_i][_qp];
}
