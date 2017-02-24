/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef BODYFORCECOMPONENT_H
#define BODYFORCECOMPONENT_H

#include "Kernel.h"

class BodyForceComponent;

template<>
InputParameters validParams<BodyForceComponent>();

class BodyForceComponent : public Kernel
{
public:
  BodyForceComponent(const InputParameters & parameters);

  virtual ~BodyForceComponent(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Coupled variables
  const VariableValue & _space_charge_density;

  // Gradients
  const VariableGradient & _grad_potential;

  // Variable numberings
  unsigned _space_charge_density_var_number;
  unsigned _potential_var_number;

  // Parameters
  unsigned _component;
  Real _p;
};

#endif // BODYFORCECOMPONENT_H
