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
#ifndef SPALARTALLMARASNOBCBC_H
#define SPALARTALLMARASNOBCBC_H

#include "IntegratedBC.h"

class SpalartAllmarasNoBCBC;

template<>
InputParameters validParams<SpalartAllmarasNoBCBC>();

class SpalartAllmarasNoBCBC : public IntegratedBC
{
public:

  SpalartAllmarasNoBCBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

private:
  Real _rho;
  Real _mu_mol;
};

#endif
