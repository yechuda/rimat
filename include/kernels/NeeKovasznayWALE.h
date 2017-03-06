/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef NEEKOVASZNAYWALE_H
#define NEEKOVASZNAYWALE_H

#include "Kernel.h"

// Forward Declarations
class NeeKovasznayWALE;

template<>
InputParameters validParams<NeeKovasznayWALE>();

class NeeKovasznayWALE : public Kernel
{
public:
  NeeKovasznayWALE(const InputParameters & parameters);

  virtual ~NeeKovasznayWALE(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Coupled variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
  const VariableValue & _d;

  // Coupled gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;

  // Variable numberings
  unsigned _u_vel_var_number;
  unsigned _v_vel_var_number;
  unsigned _w_vel_var_number;
  unsigned _d_var_number;

  // Required parameters
  Real _rho;
  Real _mu_mol;
  Real _A;
  Real _B;
};

#endif
