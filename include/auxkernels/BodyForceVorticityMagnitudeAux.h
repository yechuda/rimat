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

#ifndef BODYFORCEVORTICITYMAGNITUDEAUX_H
#define BODYFORCEVORTICITYMAGNITUDEAUX_H

#include "AuxKernel.h"

class BodyForceVorticityMagnitudeAux;

template<>
InputParameters validParams<BodyForceVorticityMagnitudeAux>();

class BodyForceVorticityMagnitudeAux : public AuxKernel
{
public:
  BodyForceVorticityMagnitudeAux(const InputParameters & parameters);

  virtual ~BodyForceVorticityMagnitudeAux() {}

protected:

  virtual Real computeValue();

  // Coupled variables
  const VariableGradient & _grad_body_force_x;
  const VariableGradient & _grad_body_force_y;
  const VariableGradient & _grad_body_force_z;
};

#endif //BODYFORCEVORTICITYMAGNITUDEAUX_H
