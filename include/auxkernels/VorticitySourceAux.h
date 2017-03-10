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

#ifndef VORTICITYSOURCEAUX_H
#define VORTICITYSOURCEAUX_H

#include "AuxKernel.h"

class VorticitySourceAux;

template<>
InputParameters validParams<VorticitySourceAux>();

class VorticitySourceAux : public AuxKernel
{
public:

  VorticitySourceAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

private:

  const VariableGradient & _grad_potential;
  const VariableGradient & _grad_density;
};

#endif // VORTICITYSOURCEAUX_H
