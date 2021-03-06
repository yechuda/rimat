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
  // Real vol = _current_elem_volume;
  // Real h = std::pow(vol, 0.33333333);
  Real h = 0.001;

  return _mu_mol + _rho * std::pow(h, 2.0) * std::abs(_grad_v_old[_qp](0));
}
