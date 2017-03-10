/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef NEEKOVASZNAYRZ_H
#define NEEKOVASZNAYRZ_H

#include "Kernel.h"

// Forward Declarations
class NeeKovasznayRZ;

template<>
InputParameters validParams<NeeKovasznayRZ>();

class NeeKovasznayRZ : public Kernel
{
public:
  NeeKovasznayRZ(const InputParameters & parameters);

  virtual ~NeeKovasznayRZ(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
};

#endif
