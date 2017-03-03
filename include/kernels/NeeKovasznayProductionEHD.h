/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef NEEKOVASZNAYPRODUCTIONEHD_H
#define NEEKOVASZNAYPRODUCTIONEHD_H

#include "Kernel.h"

// Forward Declarations
class NeeKovasznayProductionEHD;

template<>
InputParameters validParams<NeeKovasznayProductionEHD>();

class NeeKovasznayProductionEHD : public Kernel
{
public:
  NeeKovasznayProductionEHD(const InputParameters & parameters);

  virtual ~NeeKovasznayProductionEHD(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Coupled gradients
  const VariableGradient & _grad_body_force_x;
  const VariableGradient & _grad_body_force_y;
  const VariableGradient & _grad_body_force_z;

  // Variable numberings
  unsigned _body_force_x_var;
  unsigned _body_force_y_var;
  unsigned _body_force_z_var;

  // Required parameters
  Real _rho;
  Real _mu_mol;
  Real _C;
};

#endif
