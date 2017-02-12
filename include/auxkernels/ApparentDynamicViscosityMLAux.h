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

#ifndef APPARENTDYNAMICVISCOSITYMLAUX_H
#define APPARENTDYNAMICVISCOSITYMLAUX_H

#include "AuxKernel.h"

class ApparentDynamicViscosityMLAux;

template<>
InputParameters validParams<ApparentDynamicViscosityMLAux>();

class ApparentDynamicViscosityMLAux : public AuxKernel
{
public:
  ApparentDynamicViscosityMLAux(const InputParameters & parameters);

  virtual ~ApparentDynamicViscosityMLAux() {}

protected:

  virtual Real computeValue();

  // Coupled variables
  const VariableGradient & _grad_v_old;

  // Required parameters
  Real _mu_mol;
  Real _rho;
  Real _coef;
  Real _orig;
};

#endif //APPARENTDYNAMICVISCOSITYMLAUX_H
