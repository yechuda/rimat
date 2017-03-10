/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMTRACTIONFORMRANSRZ_H
#define INSMOMENTUMTRACTIONFORMRANSRZ_H

#include "INSMomentumTractionFormRANS.h"

// Forward Declarations
class INSMomentumTractionFormRANSRZ;

template<>
InputParameters validParams<INSMomentumTractionFormRANSRZ>();

/**
 * This class computes additional momentum equation residual and
 * Jacobian contributions for the incompressible Navier-Stokes
 * momentum equation in RZ (axisymmetric cylindrical) coordinates.
 */
class INSMomentumTractionFormRANSRZ : public INSMomentumTractionFormRANS
{
public:
  INSMomentumTractionFormRANSRZ(const InputParameters & parameters);

  virtual ~INSMomentumTractionFormRANSRZ(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);
};


#endif
