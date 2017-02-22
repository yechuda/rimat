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

#ifndef WALLDISTANCEAUX_H
#define WALLDISTANCEAUX_H

#include "AuxKernel.h"

class WallDistanceAux;

template<>
InputParameters validParams<WallDistanceAux>();

class WallDistanceAux : public AuxKernel
{
public:
  WallDistanceAux(const InputParameters & parameters);

  virtual ~WallDistanceAux() {}

protected:

  virtual Real computeValue();

  // Coupled variables
  const VariableValue & _G;

  // Required parameters
  Real _G0;
};

#endif //WALLDISTANCEAUX_H
