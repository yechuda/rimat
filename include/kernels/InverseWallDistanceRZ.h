/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INVERSEWALLDISTANCERZ_H
#define INVERSEWALLDISTANCERZ_H

#include "Kernel.h"

class InverseWallDistanceRZ;

template<>
InputParameters validParams<InverseWallDistanceRZ>();

class InverseWallDistanceRZ : public Kernel
{
public:
  InverseWallDistanceRZ(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  const Real _sigma;
};

#endif //INVERSEWALLDISTANCERZ_H
