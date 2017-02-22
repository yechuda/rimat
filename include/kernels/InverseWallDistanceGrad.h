/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INVERSEWALLDISTANCEGRAD_H
#define INVERSEWALLDISTANCEGRAD_H

#include "Kernel.h"

class InverseWallDistanceGrad;

template<>
InputParameters validParams<InverseWallDistanceGrad>();

class InverseWallDistanceGrad : public Kernel
{
public:
  InverseWallDistanceGrad(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:

};

#endif //INVERSEWALLDISTANCEGRAD_H
