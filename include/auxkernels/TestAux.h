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

#ifndef TESTAUX_H
#define TESTAUX_H

#include "AuxKernel.h"

class TestAux;

template<>
InputParameters validParams<TestAux>();

class TestAux : public AuxKernel
{
public:
  TestAux(const InputParameters & parameters);

  virtual ~TestAux() {}

protected:

  virtual Real computeValue();

  // Coupled variables
  const VariableGradient & _grad_u_old;
  const VariableGradient & _grad_v_old;

  // Required parameters
  Real _mu_mol;
  Real _rho;
  Real _Cs;

  const Real & _current_neighbor_volume;
};

#endif //TESTAUX_H
