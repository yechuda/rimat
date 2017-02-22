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

#ifndef APPARENTDYNAMICVISCOSITYWALEAVERAGEAUX_H
#define APPARENTDYNAMICVISCOSITYWALEAVERAGEAUX_H

#include "AuxKernel.h"

class ApparentDynamicViscosityWALEAverageAux;

template<>
InputParameters validParams<ApparentDynamicViscosityWALEAverageAux>();

class ApparentDynamicViscosityWALEAverageAux : public AuxKernel
{
public:
  ApparentDynamicViscosityWALEAverageAux(const InputParameters & parameters);

  virtual ~ApparentDynamicViscosityWALEAverageAux() {}

protected:

  virtual Real computeValue();

  // Coupled variables
  const VariableGradient & _grad_u_old;
  const VariableGradient & _grad_v_old;
  const VariableGradient & _grad_w_old;
  const VariableValue & _length_scale;

  // Required parameters
  Real _mu_mol;
  Real _rho;
  Real _Cs;

  // Old value
  const VariableValue & _mu_old;
};

#endif //APPARENTDYNAMICVISCOSITYWALEAVERAGEAUX_H
