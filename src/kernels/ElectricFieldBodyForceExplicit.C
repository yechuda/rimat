/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ElectricFieldBodyForceExplicit.h"

template<>
InputParameters validParams<ElectricFieldBodyForceExplicit>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("body_force_x", "The component of body force in x direction");
  params.addCoupledVar("body_force_y", 0, "The component of body force in y direction"); // only required in 2D and 3D
  params.addCoupledVar("body_force_z", 0, "The component of body force in z direction"); // only required in 3D

  // Required parameters
  params.addRequiredParam<unsigned>("component", "0,1,2 depending on if we are solving the x,y,z component of the momentum equation");
  params.addParam<Real>("scaling", 1.0, "Scaling factor for body force");

  return params;
}



ElectricFieldBodyForceExplicit::ElectricFieldBodyForceExplicit(const InputParameters & parameters) :
  Kernel(parameters),

  // Coupled variables
  _body_force_x(coupledValue("body_force_x")),
  _body_force_y(coupledValue("body_force_y")),
  _body_force_z(coupledValue("body_force_z")),

  // Variable numberings
  _body_force_x_var_number(coupled("body_force_x")),
  _body_force_y_var_number(coupled("body_force_y")),
  _body_force_z_var_number(coupled("body_force_z")),

  // Required parameters
  _component(getParam<unsigned>("component")),
  _scaling(getParam<Real>("scaling"))

{
}

Real ElectricFieldBodyForceExplicit::computeQpResidual()
{
  RealVectorValue body_force(_body_force_x[_qp], _body_force_y[_qp], _body_force_z[_qp]);
  return -_scaling * body_force(_component) * _test[_i][_qp];
}

Real ElectricFieldBodyForceExplicit::computeQpJacobian()
{
  return 0;
}

Real ElectricFieldBodyForceExplicit::computeQpOffDiagJacobian(unsigned jvar)
{
  if ((jvar == _body_force_x_var_number) || (jvar == _body_force_y_var_number) || (jvar == _body_force_z_var_number))
  {
    return -_scaling * _phi[_j][_qp] * _test[_i][_qp];
  }

  else
    return 0;
}
