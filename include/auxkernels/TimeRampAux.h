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

#ifndef TIMERAMPAUX_H
#define TIMERAMPAUX_H

#include "AuxKernel.h"

class TimeRampAux;

template<>
InputParameters validParams<TimeRampAux>();

class TimeRampAux : public AuxKernel
{
public:
  TimeRampAux(const InputParameters & parameters);

  virtual ~TimeRampAux() {}

protected:

  virtual Real computeValue();

  const VariableValue & _var_to_ramp;

  Real _initial_scaling;
  Real _ramp_factor;

};

#endif //TIMERAMPAUX_H
