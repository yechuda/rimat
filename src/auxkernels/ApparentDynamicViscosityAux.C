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

#include "ApparentDynamicViscosityAux.h"

template<>
InputParameters validParams<ApparentDynamicViscosityAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables
  params.addRequiredCoupledVar("nu", "apparent kinematic viscosity");

  // Required parameters
  params.addRequiredParam<Real>("rho", "density");

  return params;
}

ApparentDynamicViscosityAux::ApparentDynamicViscosityAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    // Coupled variables
    _nu(coupledValue("nu")),

    // Required parameters
    _rho(getParam<Real>("rho"))

{
}

Real ApparentDynamicViscosityAux::computeValue()
{
  return _nu[_qp] * _rho;
}
