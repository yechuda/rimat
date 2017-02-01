/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef BODYFORCEDUMMY_H
#define BODYFORCEDUMMY_H

#include "Kernel.h"

class BodyForceDummy;

template<>
InputParameters validParams<BodyForceDummy>();

class BodyForceDummy : public Kernel
{
public:
  BodyForceDummy(const InputParameters & parameters);

  virtual ~BodyForceDummy(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Parameters
  Real _value;
};

#endif // BODYFORCEDUMMY_H
