/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ELECTRICFIELDBODYFORCE_H
#define ELECTRICFIELDBODYFORCE_H

#include "Kernel.h"

class ElectricFieldBodyForce;

template<>
InputParameters validParams<ElectricFieldBodyForce>();

class ElectricFieldBodyForce : public Kernel
{
public:
  ElectricFieldBodyForce(const InputParameters & parameters);

  virtual ~ElectricFieldBodyForce(){}

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
};

#endif // ELECTRICFIELDBODYFORCE_H
