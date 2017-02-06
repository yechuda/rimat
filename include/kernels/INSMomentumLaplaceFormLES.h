/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMLAPLACEFORMLES_H
#define INSMOMENTUMLAPLACEFORMLES_H

#include "INSMomentumBaseLES.h"

// Forward Declarations
class INSMomentumLaplaceFormLES;

template<>
InputParameters validParams<INSMomentumLaplaceFormLES>();

/**
 * This class computes momentum equation residual and Jacobian viscous
 * contributions for the "Laplacian" form of the governing equations.
 */
class INSMomentumLaplaceFormLES : public INSMomentumBaseLES
{
public:
  INSMomentumLaplaceFormLES(const InputParameters & parameters);

  virtual ~INSMomentumLaplaceFormLES(){}

protected:
  virtual Real computeQpResidualViscousPart() override;
  virtual Real computeQpJacobianViscousPart() override;
  virtual Real computeQpOffDiagJacobianViscousPart(unsigned jvar) override;
};


#endif
