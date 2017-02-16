/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSEPSILONTRANSPORT_H
#define INSEPSILONTRANSPORT_H

#include "Kernel.h"

// Forward Declarations
class INSepsilonTransport;

template<>
InputParameters validParams<INSepsilonTransport>();

class INSepsilonTransport : public Kernel
{
public:
  INSepsilonTransport(const InputParameters & parameters);

  virtual ~INSepsilonTransport(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Coupled variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
  const VariableValue & _k;

  // Gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;

  // Variable numberings
  unsigned _u_vel_var_number;
  unsigned _v_vel_var_number;
  unsigned _w_vel_var_number;
  unsigned _k_var_number;

  // Parameters
  Real _Cmu;
  Real _sigma_epsilon;
  Real _C_epsilon1;
  Real _C_epsilon2;
};


#endif
