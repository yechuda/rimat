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

#ifndef APPARENTDYNAMICVISCOSITYNIKURADSEAUX_H
#define APPARENTDYNAMICVISCOSITYNIKURADSEAUX_H

#include "AuxKernel.h"

class ApparentDynamicViscosityNikuradseAux;

template<>
InputParameters validParams<ApparentDynamicViscosityNikuradseAux>();

class ApparentDynamicViscosityNikuradseAux : public AuxKernel
{
public:
  ApparentDynamicViscosityNikuradseAux(const InputParameters & parameters);

  virtual ~ApparentDynamicViscosityNikuradseAux() {}

protected:

  virtual Real computeValue();

  // Coupled variables
  const VariableGradient & _grad_v_old;

  // Required parameters
  Real _mu_mol;
  Real _rho;
  Real _D;
};

#endif //APPARENTDYNAMICVISCOSITYNIKURADSEAUX_H
