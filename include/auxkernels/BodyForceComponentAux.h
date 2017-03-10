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

#ifndef BODYFORCECOMPONENTAUX_H
#define BODYFORCECOMPONENTAUX_H

#include "AuxKernel.h"

class BodyForceComponentAux;

template<>
InputParameters validParams<BodyForceComponentAux>();

class BodyForceComponentAux : public AuxKernel
{
public:

  BodyForceComponentAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

private:

  const VariableGradient & _grad_potential;
  const VariableValue & _space_charge_density;
  int _component;
};

#endif // BODYFORCECOMPONENTAUX_H
