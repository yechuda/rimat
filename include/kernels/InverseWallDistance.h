/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INVERSEWALLDISTANCE_H
#define INVERSEWALLDISTANCE_H

#include "Kernel.h"

class InverseWallDistance;

template<>
InputParameters validParams<InverseWallDistance>();

class InverseWallDistance : public Kernel
{
public:
  InverseWallDistance(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const Real _sigma;
};

#endif //INVERSEWALLDISTANCE_H
