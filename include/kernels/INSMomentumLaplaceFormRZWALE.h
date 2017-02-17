/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMLAPLACEFORMRZWALE_H
#define INSMOMENTUMLAPLACEFORMRZWALE_H

#include "INSMomentumLaplaceFormWALE.h"

// Forward Declarations
class INSMomentumLaplaceFormRZWALE;

template<>
InputParameters validParams<INSMomentumLaplaceFormRZWALE>();

/**
 * This class computes additional momentum equation residual and
 * Jacobian contributions for the incompressible Navier-Stokes
 * momentum equation in RZ (axisymmetric cylindrical) coordinates,
 * using the "Laplace" form of the governing equations.
 */
class INSMomentumLaplaceFormRZWALE : public INSMomentumLaplaceFormWALE
{
public:
  INSMomentumLaplaceFormRZWALE(const InputParameters & parameters);

  virtual ~INSMomentumLaplaceFormRZWALE(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);
};


#endif
