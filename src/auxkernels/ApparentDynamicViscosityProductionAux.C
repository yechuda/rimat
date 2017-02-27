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

#include "ApparentDynamicViscosityProductionAux.h"

template<>
InputParameters validParams<ApparentDynamicViscosityProductionAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables
  params.addRequiredCoupledVar("d", "wall distance");
  params.addRequiredCoupledVar("body_force_vorticity_mag", "The vorticity of the body force");

  // Required parameters
  params.addRequiredParam<Real>("mu_mol", "molecular dynamic viscosity");
  params.addRequiredParam<Real>("rho", "density");
  params.addRequiredParam<Real>("Ce", "scalling coefficient");

  return params;
}

ApparentDynamicViscosityProductionAux::ApparentDynamicViscosityProductionAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    // Coupled variables
    _d(coupledValue("d")),
    _body_force_vorticity_mag(coupledValue("body_force_vorticity_mag")),

    // Required parameters
    _mu_mol(getParam<Real>("mu_mol")),
    _rho(getParam<Real>("rho")),
    _Ce(getParam<Real>("Ce"))

{
}

Real ApparentDynamicViscosityProductionAux::computeValue()
{
  return _mu_mol + _Ce * std::pow(_body_force_vorticity_mag[_qp] / _rho, 0.5) * std::pow(_d[_qp], 2.0) * _rho;
}
