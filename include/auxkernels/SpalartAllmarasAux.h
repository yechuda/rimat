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

#ifndef SPALARTALLMARASAUX_H
#define SPALARTALLMARASAUX_H

#include "AuxKernel.h"

class SpalartAllmarasAux;

template<>
InputParameters validParams<SpalartAllmarasAux>();

class SpalartAllmarasAux : public AuxKernel
{
public:
  SpalartAllmarasAux(const InputParameters & parameters);

  virtual ~SpalartAllmarasAux() {}

protected:

  virtual Real computeValue();

  // Coupled variables
  const VariableValue & _nu_tilde;

  // Required parameters
  Real _mu_mol;
  Real _rho;
};

#endif //SPALARTALLMARASAUX_H
