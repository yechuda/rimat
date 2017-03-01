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
#ifndef NEEKOVASZNAYNOBCBC_H
#define NEEKOVASZNAYNOBCBC_H

#include "IntegratedBC.h"

class NeeKovasznayNoBCBC;

template<>
InputParameters validParams<NeeKovasznayNoBCBC>();

class NeeKovasznayNoBCBC : public IntegratedBC
{
public:

  NeeKovasznayNoBCBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

};

#endif
