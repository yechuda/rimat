/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef NEEKOVASZNAYPRODUCTIONEHDSIMPLE_H
#define NEEKOVASZNAYPRODUCTIONEHDSIMPLE_H

#include "Kernel.h"

// Forward Declarations
class NeeKovasznayProductionEHDSimple;

template<>
InputParameters validParams<NeeKovasznayProductionEHDSimple>();

class NeeKovasznayProductionEHDSimple : public Kernel
{
public:
  NeeKovasznayProductionEHDSimple(const InputParameters & parameters);

  virtual ~NeeKovasznayProductionEHDSimple(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Coupled variables
  const VariableValue & _d;

  // Coupled gradients
  const VariableGradient & _grad_body_force_x;
  const VariableGradient & _grad_body_force_y;
  const VariableGradient & _grad_body_force_z;

  // Variable numberings
  unsigned _body_force_x_var;
  unsigned _body_force_y_var;
  unsigned _body_force_z_var;
  unsigned _d_var;

  // Required parameters
  Real _rho;
  Real _C;
};

#endif
