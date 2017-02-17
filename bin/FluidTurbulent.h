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

#ifndef FLUIDTURBULENT_H
#define FLUIDTURBULENT_H

#include "Material.h"

//Forward Declarations
class FluidTurbulent;

template<>
InputParameters validParams<FluidTurbulent>();

class FluidTurbulent : public Material
{
public:
  FluidTurbulent(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

private:

  // Declarations
  MaterialProperty<Real> & _rho;
  MaterialProperty<Real> & _mu;

  // Gradients
  const VariableGradient & _grad_u;
  const VariableGradient & _grad_v;
  const VariableGradient & _grad_w;

  // Required parameters
  Real _mu_mol;
  Real _rho_param;
  Real _Cs;

};

#endif //FLUIDTURBULENT_H
