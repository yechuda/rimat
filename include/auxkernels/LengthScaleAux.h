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

#ifndef LENGTHSCALEAUX_H
#define LENGTHSCALEAUX_H

#include "AuxKernel.h"

class LengthScaleAux;

template<>
InputParameters validParams<LengthScaleAux>();

class LengthScaleAux : public AuxKernel
{
public:
  LengthScaleAux(const InputParameters & parameters);

  virtual ~LengthScaleAux() {}

protected:

  virtual Real computeValue();

  Real _D;
};

#endif //LENGTHSCALEAUX_H
