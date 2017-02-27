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

#include "TimeRampAux.h"

template<>
InputParameters validParams<TimeRampAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("var_to_ramp", "variable to ramp");
  params.addParam<Real>("initial_scaling", 1.0e-06, "Initial value of scaling factor");
  params.addParam<Real>("ramp_factor", 0.04, "Multiplication factor for ramp pattern");

  return params;
}

TimeRampAux::TimeRampAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    _var_to_ramp(coupledValue("var_to_ramp")),
    _initial_scaling(getParam<Real>("initial_scaling")),
    _ramp_factor(getParam<Real>("ramp_factor"))
{
}

Real TimeRampAux::computeValue()
{
  Real _scaling = _initial_scaling + (_t_step - 1) * _ramp_factor;

  if (_scaling > 1.0)
    _scaling = 1.0;

  return _scaling * _var_to_ramp[_qp];
}
