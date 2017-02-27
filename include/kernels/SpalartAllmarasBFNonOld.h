/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef SPALARTALLMARASBFNONOLD_H
#define SPALARTALLMARASBFNONOLD_H

#include "Kernel.h"

// Forward Declarations
class SpalartAllmarasBFNonOld;

template<>
InputParameters validParams<SpalartAllmarasBFNonOld>();

class SpalartAllmarasBFNonOld : public Kernel
{
public:
  SpalartAllmarasBFNonOld(const InputParameters & parameters);

  virtual ~SpalartAllmarasBFNonOld(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  // Coupled old variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
  const VariableValue & _d;
  const VariableValue & _body_force_vorticity_mag;

  // Coupled old gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;

  // Required parameters
  Real _rho;
  Real _mu_mol;
  Real _Ce;

  // Old value
  const VariableValue & _nu_tilde;
};

#endif
