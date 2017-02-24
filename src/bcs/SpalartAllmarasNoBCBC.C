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
#include "SpalartAllmarasNoBCBC.h"

template<>
InputParameters validParams<SpalartAllmarasNoBCBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("rho", "density");
  params.addRequiredParam<Real>("mu_mol", "molecular dynamic viscosiyty");

  return params;
}

SpalartAllmarasNoBCBC::SpalartAllmarasNoBCBC(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _rho(getParam<Real>("rho")),
    _mu_mol(getParam<Real>("mu_mol"))
{}

Real
SpalartAllmarasNoBCBC::computeQpResidual()
{
  Real nu_mol = _mu_mol / _rho;
  Real sigma = 2.0 / 3.0;

  Real first_part = -nu_mol / sigma * (_grad_u[_qp] * _normals[_qp]) * _test[_i][_qp];
  Real second_part = -1.0 / sigma * _u[_qp] * (_grad_u[_qp] * _normals[_qp]) * _test[_i][_qp];
  return first_part + second_part;
}

Real
SpalartAllmarasNoBCBC::computeQpJacobian()
{
  Real nu_mol = _mu_mol / _rho;
  Real sigma = 2.0 / 3.0;

  Real first_part = -nu_mol / sigma * (_grad_phi[_j][_qp] * _normals[_qp]) * _test[_i][_qp];
  Real second_part = -1.0 / sigma * _u[_qp] * (_grad_phi[_j][_qp] * _normals[_qp]) * _test[_i][_qp] -
                      1.0 / sigma * _phi[_j][_qp] * (_grad_u[_qp] * _normals[_qp]) * _test[_i][_qp];
  return first_part + second_part;
}
