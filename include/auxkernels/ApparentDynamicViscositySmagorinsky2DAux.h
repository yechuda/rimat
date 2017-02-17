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

#ifndef APPARENTDYNAMICVISCOSITYSMAGORINSKY2DAUX_H
#define APPARENTDYNAMICVISCOSITYSMAGORINSKY2DAUX_H

#include "AuxKernel.h"

class ApparentDynamicViscositySmagorinsky2DAux;

template<>
InputParameters validParams<ApparentDynamicViscositySmagorinsky2DAux>();

class ApparentDynamicViscositySmagorinsky2DAux : public AuxKernel
{
public:
  ApparentDynamicViscositySmagorinsky2DAux(const InputParameters & parameters);

  virtual ~ApparentDynamicViscositySmagorinsky2DAux() {}

protected:

  virtual Real computeValue();

  // Coupled variables
  const VariableGradient & _grad_u_old;
  const VariableGradient & _grad_v_old;

  // Required parameters
  Real _mu_mol;
  Real _rho;
  Real _Cs;
};

#endif //APPARENTDYNAMICVISCOSITYSMAGORINSKY2DAUX_H
