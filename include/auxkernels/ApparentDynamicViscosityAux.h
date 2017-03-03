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

#ifndef APPARENTDYNAMICVISCOSITYAUX_H
#define APPARENTDYNAMICVISCOSITYAUX_H

#include "AuxKernel.h"

class ApparentDynamicViscosityAux;

template<>
InputParameters validParams<ApparentDynamicViscosityAux>();

class ApparentDynamicViscosityAux : public AuxKernel
{
public:
  ApparentDynamicViscosityAux(const InputParameters & parameters);

  virtual ~ApparentDynamicViscosityAux() {}

protected:

  virtual Real computeValue();

  // Coupled variables
  const VariableValue & _nu;

  // Required parameters
  Real _rho;
};

#endif //APPARENTDYNAMICVISCOSITYAUX_H
