/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
#ifndef BODYFORCEBC_H
#define BODYFORCEBC_H

#include "IntegratedBC.h"

class BodyForceBC;
class Function;

template<>
InputParameters validParams<BodyForceBC>();

class BodyForceBC : public IntegratedBC
{
public:

  BodyForceBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

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

#endif
