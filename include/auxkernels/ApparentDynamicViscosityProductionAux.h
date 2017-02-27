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

#ifndef APPARENTDYNAMICVISCOSITYPRODUCTIONAUX_H
#define APPARENTDYNAMICVISCOSITYPRODUCTIONAUX_H

#include "AuxKernel.h"

class ApparentDynamicViscosityProductionAux;

template<>
InputParameters validParams<ApparentDynamicViscosityProductionAux>();

class ApparentDynamicViscosityProductionAux : public AuxKernel
{
public:
  ApparentDynamicViscosityProductionAux(const InputParameters & parameters);

  virtual ~ApparentDynamicViscosityProductionAux() {}

protected:

  virtual Real computeValue();

  // Coupled variables
  const VariableValue & _d;
  const VariableValue & _body_force_vorticity_mag;

  // Required parameters
  Real _mu_mol;
  Real _rho;
  Real _Ce;
};

#endif //APPARENTDYNAMICVISCOSITYPRODUCTIONAUX_H
