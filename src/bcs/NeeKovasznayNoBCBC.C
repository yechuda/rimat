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
#include "NeeKovasznayNoBCBC.h"

template<>
InputParameters validParams<NeeKovasznayNoBCBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  return params;
}

NeeKovasznayNoBCBC::NeeKovasznayNoBCBC(const InputParameters & parameters) :
    IntegratedBC(parameters)

{}

Real
NeeKovasznayNoBCBC::computeQpResidual()
{
  return -_u[_qp] * (_grad_u[_qp] * _normals[_qp]) * _test[_i][_qp];
}

Real
NeeKovasznayNoBCBC::computeQpJacobian()
{
  return -_phi[_j][_qp] * (_grad_u[_qp] * _normals[_qp]) * _test[_i][_qp] -
          _u[_qp] * (_grad_phi[_j][_qp] * _normals[_qp]) * _test[_i][_qp];
}
