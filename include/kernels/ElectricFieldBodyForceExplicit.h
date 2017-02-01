/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ELECTRICFIELDBODYFORCEEXPLICIT_H
#define ELECTRICFIELDBODYFORCEEXPLICIT_H

#include "Kernel.h"

class ElectricFieldBodyForceExplicit;

template<>
InputParameters validParams<ElectricFieldBodyForceExplicit>();

class ElectricFieldBodyForceExplicit : public Kernel
{
public:
  ElectricFieldBodyForceExplicit(const InputParameters & parameters);

  virtual ~ElectricFieldBodyForceExplicit(){}

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
  Real _scaling;
};

#endif // ELECTRICFIELDBODYFORCEEXPLICIT_H
