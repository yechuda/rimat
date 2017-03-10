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

#ifndef SOURCEWALEAUX_H
#define SOURCEWALEAUX_H

#include "AuxKernel.h"

class SourceWALEAux;

template<>
InputParameters validParams<SourceWALEAux>();

class SourceWALEAux : public AuxKernel
{
public:
  SourceWALEAux(const InputParameters & parameters);

  virtual ~SourceWALEAux() {}

protected:

  virtual Real computeValue();

  // Coupled variables
  const VariableValue & _d;

  // Coupled gradients
  const VariableGradient & _grad_body_force_x;
  const VariableGradient & _grad_body_force_y;
  const VariableGradient & _grad_body_force_z;

  // Required parameters
  Real _C;
};

#endif //SOURCEWALEAUX_H
