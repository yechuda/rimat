/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMTRACTIONFORMRANS_H
#define INSMOMENTUMTRACTIONFORMRANS_H

#include "INSMomentumBaseRANS.h"

// Forward Declarations
class INSMomentumTractionFormRANS;

template<>
InputParameters validParams<INSMomentumTractionFormRANS>();

/**
 * This class computes momentum equation residual and Jacobian viscous
 * contributions for the "traction" form of the governing equations.
 */
class INSMomentumTractionFormRANS : public INSMomentumBaseRANS
{
public:
  INSMomentumTractionFormRANS(const InputParameters & parameters);

  virtual ~INSMomentumTractionFormRANS(){}

protected:
  virtual Real computeQpResidualViscousPart() override;
  virtual Real computeQpJacobianViscousPart() override;
  virtual Real computeQpOffDiagJacobianViscousPart(unsigned jvar) override;
};


#endif
