/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMLAPLACEFORMRZTURBULENT_H
#define INSMOMENTUMLAPLACEFORMRZTURBULENT_H

#include "INSMomentumLaplaceFormTurbulent.h"

// Forward Declarations
class INSMomentumLaplaceFormRZTurbulent;

template<>
InputParameters validParams<INSMomentumLaplaceFormRZTurbulent>();

/**
 * This class computes additional momentum equation residual and
 * Jacobian contributions for the incompressible Navier-Stokes
 * momentum equation in RZ (axisymmetric cylindrical) coordinates,
 * using the "Laplace" form of the governing equations.
 */
class INSMomentumLaplaceFormRZTurbulent : public INSMomentumLaplaceFormTurbulent
{
public:
  INSMomentumLaplaceFormRZTurbulent(const InputParameters & parameters);

  virtual ~INSMomentumLaplaceFormRZTurbulent(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);
};


#endif
