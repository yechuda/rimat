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

#ifndef APPARENTDYNAMICVISCOSITYWALEAUX_H
#define APPARENTDYNAMICVISCOSITYWALEAUX_H

#include "AuxKernel.h"

class ApparentDynamicViscosityWALEAux;

template<>
InputParameters validParams<ApparentDynamicViscosityWALEAux>();

class ApparentDynamicViscosityWALEAux : public AuxKernel
{
public:
  ApparentDynamicViscosityWALEAux(const InputParameters & parameters);

  virtual ~ApparentDynamicViscosityWALEAux() {}

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
};

#endif //APPARENTDYNAMICVISCOSITYWALEAUX_H
