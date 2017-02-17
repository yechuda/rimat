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

#include "ApparentDynamicViscositySmagorinsky2DAux.h"

template<>
InputParameters validParams<ApparentDynamicViscositySmagorinsky2DAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addRequiredCoupledVar("v", "y-velocity");

  // Required parameters
  params.addRequiredParam<Real>("mu_mol", "molecular dynamic viscosity");
  params.addRequiredParam<Real>("rho", "density");
  params.addRequiredParam<Real>("Cs", "Smagorinsky coefficient");

  return params;
}

ApparentDynamicViscositySmagorinsky2DAux::ApparentDynamicViscositySmagorinsky2DAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    // Coupled variables
    _grad_u_old(coupledGradientOld("u")),
    _grad_v_old(coupledGradientOld("v")),

    // Required parameters
    _mu_mol(getParam<Real>("mu_mol")),
    _rho(getParam<Real>("rho")),
    _Cs(getParam<Real>("Cs"))

{
}

Real ApparentDynamicViscositySmagorinsky2DAux::computeValue()
{
  // Real vol = _current_elem_volume;
  // Real h = std::pow(vol, 0.33333333);
  // Real h = 0.001;
  // return _mu_mol + _rho * std::pow(h, 2.0) * std::abs(_grad_v_old[_qp](0));

  Real OP_squared = 2.0 * std::pow(_grad_u_old[_qp](0), 2.0) + 2.0 * std::pow(_grad_v_old[_qp](1), 2.0) + std::pow(_grad_u_old[_qp](1) + _grad_v_old[_qp](0), 2.0);
  Real OP = std::pow(OP, 0.5);

  // Real lm = _Cs * 2.0 * _current_elem->hmax();
  Real vol = _current_elem_volume;
  Real lm = std::pow(vol, 0.33333333);

  return _mu_mol + _rho * std::pow(lm, 2.0) * OP;
}
