/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMLAPLACEFORMRZLES_H
#define INSMOMENTUMLAPLACEFORMRZLES_H

#include "INSMomentumLaplaceFormLES.h"

// Forward Declarations
class INSMomentumLaplaceFormRZLES;

template<>
InputParameters validParams<INSMomentumLaplaceFormRZLES>();

/**
 * This class computes additional momentum equation residual and
 * Jacobian contributions for the incompressible Navier-Stokes
 * momentum equation in RZ (axisymmetric cylindrical) coordinates,
 * using the "Laplace" form of the governing equations.
 */
class INSMomentumLaplaceFormRZLES : public INSMomentumLaplaceFormLES
{
public:
  INSMomentumLaplaceFormRZLES(const InputParameters & parameters);

  virtual ~INSMomentumLaplaceFormRZLES(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);
};


#endif
