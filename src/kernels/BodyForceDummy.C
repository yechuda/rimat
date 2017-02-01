/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "BodyForceDummy.h"

template<>
InputParameters validParams<BodyForceDummy>()
{
  InputParameters params = validParams<Kernel>();

  // Required parameters
  params.addRequiredParam<Real>("value", "The value of the body force");

  return params;
}

BodyForceDummy::BodyForceDummy(const InputParameters & parameters) :
  Kernel(parameters),

  // Required parameters
  _value(getParam<Real>("value"))

{
}

Real BodyForceDummy::computeQpResidual()
{
  return -_value * _test[_i][_qp];
}

Real BodyForceDummy::computeQpJacobian()
{
  return 0;
}

Real BodyForceDummy::computeQpOffDiagJacobian(unsigned jvar)
{
  return 0;
}
