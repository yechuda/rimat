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
#include "BodyForceComponentAux.h"

template<>
InputParameters validParams<BodyForceComponentAux>()
{
  MooseEnum component("x=0 y=1 z=2");
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("potential", "The coupled variable of potential");
  params.addRequiredCoupledVar("space_charge_density", "The coupled variable of space charge density");
  params.addParam<MooseEnum>("component", component, "The component of the body force to compute");
  return params;
}

BodyForceComponentAux::BodyForceComponentAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    _grad_potential(coupledGradient("potential")),
    _space_charge_density(coupledValue("space_charge_density")),
    _component(getParam<MooseEnum>("component"))
{
}

Real
BodyForceComponentAux::computeValue()
{
  return -_space_charge_density[_qp] * _grad_potential[_qp](_component);
}
