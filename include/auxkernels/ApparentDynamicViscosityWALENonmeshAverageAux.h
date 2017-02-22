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

#ifndef APPARENTDYNAMICVISCOSITYWALENONMESHAVERAGEAUX_H
#define APPARENTDYNAMICVISCOSITYWALENONMESHAVERAGEAUX_H

#include "AuxKernel.h"

class ApparentDynamicViscosityWALENonmeshAverageAux;

template<>
InputParameters validParams<ApparentDynamicViscosityWALENonmeshAverageAux>();

class ApparentDynamicViscosityWALENonmeshAverageAux : public AuxKernel
{
public:
  ApparentDynamicViscosityWALENonmeshAverageAux(const InputParameters & parameters);

  virtual ~ApparentDynamicViscosityWALENonmeshAverageAux() {}

protected:

  virtual Real computeValue();

  // Coupled variables
  const VariableGradient & _grad_u_old;
  const VariableGradient & _grad_v_old;
  const VariableGradient & _grad_w_old;

  // Required parameters
  Real _mu_mol;
  Real _rho;
  Real _Cs;
  Real _delta;

  // Old value
  const VariableValue & _mu_old;
};

#endif //APPARENTDYNAMICVISCOSITYWALENONMESHAVERAGEAUX_H
