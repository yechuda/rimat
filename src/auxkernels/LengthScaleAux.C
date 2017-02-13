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

#include "LengthScaleAux.h"

template<>
InputParameters validParams<LengthScaleAux>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredParam<Real>("D", "pipe density");

  return params;
}

LengthScaleAux::LengthScaleAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    _D(getParam<Real>("D"))

{
}

Real LengthScaleAux::computeValue()
{
  Real r = _q_point[_qp](0);
  Real R = _D / 2.0;

  if (r > R)
    r = R;

  Real lm_squared;

  lm_squared = std::pow(R, 2.0) * (0.03125 - 0.03125 * std::pow(r / R, 2.0));

  return lm_squared;
}
