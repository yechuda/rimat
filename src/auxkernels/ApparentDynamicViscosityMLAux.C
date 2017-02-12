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

#include "ApparentDynamicViscosityMLAux.h"

template<>
InputParameters validParams<ApparentDynamicViscosityMLAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables
  params.addRequiredCoupledVar("v", "y-velocity");

  // Required parameters
  params.addRequiredParam<Real>("mu_mol", "molecular dynamic viscosity");
  params.addRequiredParam<Real>("rho", "density");
  params.addRequiredParam<Real>("coef", "coefficient");
  params.addRequiredParam<Real>("jet_origin_y", "y coordinate of jet origin");

  return params;
}

ApparentDynamicViscosityMLAux::ApparentDynamicViscosityMLAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    // Coupled variables
    _grad_v_old(coupledGradientOld("v")),

    // Required parameters
    _mu_mol(getParam<Real>("mu_mol")),
    _rho(getParam<Real>("rho")),
    _coef(getParam<Real>("coef")),
    _orig(getParam<Real>("jet_origin_y"))

{
}

Real ApparentDynamicViscosityMLAux::computeValue()
{
  Real y = _q_point[_qp](1);
  Real pos;

  if (y < _orig)
    {
      pos = 0;
    }
  else
    {
      pos = y - _orig;
    }

  Real ml = _coef * pos;

  return _mu_mol + _rho * std::pow(ml, 2.0) * std::abs(_grad_v_old[_qp](0));
}
