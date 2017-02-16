/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMLAPLACEFORMWALE_H
#define INSMOMENTUMLAPLACEFORMWALE_H

#include "INSMomentumBaseWALE.h"

// Forward Declarations
class INSMomentumLaplaceFormWALE;

template<>
InputParameters validParams<INSMomentumLaplaceFormWALE>();

/**
 * This class computes momentum equation residual and Jacobian viscous
 * contributions for the "Laplacian" form of the governing equations.
 */
class INSMomentumLaplaceFormWALE : public INSMomentumBaseWALE
{
public:
  INSMomentumLaplaceFormWALE(const InputParameters & parameters);

  virtual ~INSMomentumLaplaceFormWALE(){}

protected:
  virtual Real computeQpResidualViscousPart() override;
  virtual Real computeQpJacobianViscousPart() override;
  virtual Real computeQpOffDiagJacobianViscousPart(unsigned jvar) override;
};


#endif
