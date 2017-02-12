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

#include "ApparentDynamicViscosityNikuradseAux.h"

template<>
InputParameters validParams<ApparentDynamicViscosityNikuradseAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables
  params.addRequiredCoupledVar("v", "y-velocity");

  // Required parameters
  params.addRequiredParam<Real>("mu_mol", "molecular dynamic viscosity");
  params.addRequiredParam<Real>("rho", "density");
  params.addRequiredParam<Real>("D", "pipe diameter");

  return params;
}

ApparentDynamicViscosityNikuradseAux::ApparentDynamicViscosityNikuradseAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    // Coupled variables
    _grad_v_old(coupledGradientOld("v")),

    // Required parameters
    _mu_mol(getParam<Real>("mu_mol")),
    _rho(getParam<Real>("rho")),
    _D(getParam<Real>("D"))

{
}

Real ApparentDynamicViscosityNikuradseAux::computeValue()
{
  Real r = _q_point[_qp](0);
  Real wall_dist;

  if (r > (_D / 2.0))
    {
      wall_dist = 0.0;
    }
  else
    {
      wall_dist = (_D / 2.0) - r;
    }

  Real rel_dist = 1.0 - (2.0 * wall_dist / _D);
  Real ml = (_D / 2.0) * (0.14 - 0.08 * std::pow(rel_dist, 2.0) - 0.06 * std::pow(rel_dist, 4.0));

  return _mu_mol + _rho * std::pow(ml, 2.0) * std::abs(_grad_v_old[_qp](0));
}
