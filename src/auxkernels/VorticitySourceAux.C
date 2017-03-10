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
#include "VorticitySourceAux.h"

template<>
InputParameters validParams<VorticitySourceAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("potential", "The coupled variable of potential");
  params.addRequiredCoupledVar("space_charge_density", "The coupled variable of space charge density");
  return params;
}

VorticitySourceAux::VorticitySourceAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    _grad_potential(coupledGradient("potential")),
    _grad_density(coupledGradient("space_charge_density"))
{
}

Real
VorticitySourceAux::computeValue()
{
  Real S_mag_squared = std::pow(_grad_density[_qp](2) * _grad_potential[_qp](1) - _grad_density[_qp](1) * _grad_potential[_qp](2), 2.0) +
                       std::pow(_grad_density[_qp](0) * _grad_potential[_qp](2) - _grad_density[_qp](2) * _grad_potential[_qp](0), 2.0) +
                       std::pow(_grad_density[_qp](1) * _grad_potential[_qp](0) - _grad_density[_qp](0) * _grad_potential[_qp](1), 2.0);
  Real S_mag = std::pow(S_mag_squared, 0.5);
  return S_mag;
}
