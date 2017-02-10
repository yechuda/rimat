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

#ifndef ELEMVOLUMEAUX_H
#define ELEMVOLUMEAUX_H

#include "AuxKernel.h"

class ElemVolumeAux;

template<>
InputParameters validParams<ElemVolumeAux>();

class ElemVolumeAux : public AuxKernel
{
public:
  ElemVolumeAux(const InputParameters & parameters);

  virtual ~ElemVolumeAux() {}

protected:

  virtual Real computeValue();
};

#endif //ELEMVOLUMEAUX_H
