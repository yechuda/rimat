/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ElectricFieldBodyForce.h"

template<>
InputParameters validParams<ElectricFieldBodyForce>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("potential", "electric potential");
  params.addRequiredCoupledVar("space_charge_density", "space charge density");

  // Required parameters
  params.addRequiredParam<unsigned>("component", "0,1,2 depending on if we are solving the x,y,z component of the momentum equation");

  return params;
}



ElectricFieldBodyForce::ElectricFieldBodyForce(const InputParameters & parameters) :
  Kernel(parameters),

  // Coupled variables
  _space_charge_density(coupledValue("space_charge_density")),

  // Gradients
  _grad_potential(coupledGradient("potential")),

  // Variable numberings
  _space_charge_density_var_number(coupled("space_charge_density")),
  _potential_var_number(coupled("potential")),

  // Required parameters
  _component(getParam<unsigned>("component"))

{
}

Real ElectricFieldBodyForce::computeQpResidual()
{
  return _space_charge_density[_qp] * _grad_potential[_qp](_component) * _test[_i][_qp];
}

Real ElectricFieldBodyForce::computeQpJacobian()
{
  return 0;
}

Real ElectricFieldBodyForce::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _space_charge_density_var_number)
  {
    return _phi[_j][_qp] * _grad_potential[_qp](_component) * _test[_i][_qp];
  }

  else if (jvar == _potential_var_number)
  {
    return _space_charge_density[_qp] * _grad_phi[_j][_qp](_component) * _test[_i][_qp];
  }

  else
    return 0;
}
