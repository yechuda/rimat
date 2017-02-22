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

#include "WallDistanceAux.h"

template<>
InputParameters validParams<WallDistanceAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables
  params.addRequiredCoupledVar("G", "inverse wall distance");

  // Required parameters
  params.addRequiredParam<Real>("G0", "inverse offset distance");

  return params;
}

WallDistanceAux::WallDistanceAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    // Coupled variables
    _G(coupledValue("G")),

    // Required parameters
    _G0(getParam<Real>("G0"))

{
}

Real WallDistanceAux::computeValue()
{
  return (1.0 / _G[_qp]) - (1.0 / _G0);
}
