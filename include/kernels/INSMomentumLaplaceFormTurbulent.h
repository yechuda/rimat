/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMLAPLACEFORMTURBULENT_H
#define INSMOMENTUMLAPLACEFORMTURBULENT_H

#include "INSMomentumBaseTurbulent.h"

// Forward Declarations
class INSMomentumLaplaceFormTurbulent;

template<>
InputParameters validParams<INSMomentumLaplaceFormTurbulent>();

/**
 * This class computes momentum equation residual and Jacobian viscous
 * contributions for the "Laplacian" form of the governing equations.
 */
class INSMomentumLaplaceFormTurbulent : public INSMomentumBaseTurbulent
{
public:
  INSMomentumLaplaceFormTurbulent(const InputParameters & parameters);

  virtual ~INSMomentumLaplaceFormTurbulent(){}

protected:
  virtual Real computeQpResidualViscousPart() override;
  virtual Real computeQpJacobianViscousPart() override;
  virtual Real computeQpOffDiagJacobianViscousPart(unsigned jvar) override;
};


#endif
