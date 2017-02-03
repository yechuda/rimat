/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ELECTRICFIELDBODYFORCEEXPLICITRAMP_H
#define ELECTRICFIELDBODYFORCEEXPLICITRAMP_H

#include "Kernel.h"

class ElectricFieldBodyForceExplicitRamp;

template<>
InputParameters validParams<ElectricFieldBodyForceExplicitRamp>();

class ElectricFieldBodyForceExplicitRamp : public Kernel
{
public:
  ElectricFieldBodyForceExplicitRamp(const InputParameters & parameters);

  virtual ~ElectricFieldBodyForceExplicitRamp(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Coupled variables
  const VariableValue & _body_force_x;
  const VariableValue & _body_force_y;
  const VariableValue & _body_force_z;

  // Variable numberings
  unsigned _body_force_x_var_number;
  unsigned _body_force_y_var_number;
  unsigned _body_force_z_var_number;

  // Parameters
  unsigned _component;
  Real _initial_scaling;
  Real _ramp_factor;
};

#endif // ELECTRICFIELDBODYFORCEEXPLICITRAMP_H
