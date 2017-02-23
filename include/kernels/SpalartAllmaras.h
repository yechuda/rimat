/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef SPALARTALLMARAS_H
#define SPALARTALLMARAS_H

#include "Kernel.h"

// Forward Declarations
class SpalartAllmaras;

template<>
InputParameters validParams<SpalartAllmaras>();

class SpalartAllmaras : public Kernel
{
public:
  SpalartAllmaras(const InputParameters & parameters);

  virtual ~SpalartAllmaras(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  // Coupled old variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
  const VariableValue & _d;

  // Coupled old gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;

  // Required parameters
  Real _rho;
  Real _mu_mol;

  // Old value
  const VariableValue & _nu_tilde;
};

#endif
