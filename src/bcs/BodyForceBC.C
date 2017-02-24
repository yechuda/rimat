/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
#include "BodyForceBC.h"
#include "Function.h"

template<>
InputParameters validParams<BodyForceBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredCoupledVar("potential", "electric potential");
  params.addRequiredCoupledVar("space_charge_density", "space charge density");
  params.addRequiredParam<unsigned>("component", "0,1,2 depending on if we are solving the x,y,z component of the momentum equation");
  params.addRequiredParam<Real>("penalty", "Penalty scalar");

  return params;
}

BodyForceBC::BodyForceBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    // Coupled variables
    _space_charge_density(coupledValue("space_charge_density")),

    // Gradients
    _grad_potential(coupledGradient("potential")),

    // Variable numberings
    _space_charge_density_var_number(coupled("space_charge_density")),
    _potential_var_number(coupled("potential")),

    // Required parameters
    _component(getParam<unsigned>("component")),
    _p(getParam<Real>("penalty"))
{}

Real
BodyForceBC::computeQpResidual()
{
  return _p * _test[_i][_qp] * (_space_charge_density[_qp] * _grad_potential[_qp](_component) + _u[_qp]);
}

Real
BodyForceBC::computeQpJacobian()
{
  return _p * _test[_i][_qp] * _phi[_j][_qp];
}

Real
BodyForceBC::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _space_charge_density_var_number)
  {
    return _p * _phi[_j][_qp] * _grad_potential[_qp](_component) * _test[_i][_qp];
  }

  else if (jvar == _potential_var_number)
  {
    return _p * _space_charge_density[_qp] * _grad_phi[_j][_qp](_component) * _test[_i][_qp];
  }

  else
    return 0;
}
