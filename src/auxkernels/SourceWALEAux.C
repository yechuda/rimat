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

#include "SourceWALEAux.h"

template<>
InputParameters validParams<SourceWALEAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables
  params.addRequiredCoupledVar("body_force_x", "The component of body force in x direction");
  params.addCoupledVar("body_force_y", 0, "The component of body force in y direction"); // only required in 2D and 3D
  params.addCoupledVar("body_force_z", 0, "The component of body force in z direction"); // only required in 3D
  params.addRequiredCoupledVar("d", "wall distance");

  // Required parameters
  params.addRequiredParam<Real>("C", "EHD production coefficient");

  return params;
}

SourceWALEAux::SourceWALEAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    // Coupled variables
    _d(coupledValue("d")),

    // Coupled gradients
    _grad_body_force_x(coupledGradient("body_force_x")),
    _grad_body_force_y(coupledGradient("body_force_y")),
    _grad_body_force_z(coupledGradient("body_force_z")),

    // Required parameters
    _C(getParam<Real>("C"))

{
}

Real SourceWALEAux::computeValue()
{
  RealTensorValue Sij;
  Sij(0,0) =            _grad_body_force_x[_qp](0);
  Sij(0,1) = Sij(1,0) = 0.5 * (_grad_body_force_x[_qp](1) + _grad_body_force_y[_qp](0));
  Sij(1,1) =            _grad_body_force_y[_qp](1);
  Sij(0,2) = Sij(2,0) = 0.5 * (_grad_body_force_x[_qp](2) + _grad_body_force_z[_qp](0));
  Sij(1,2) = Sij(2,1) = 0.5 * (_grad_body_force_y[_qp](2) + _grad_body_force_z[_qp](1));
  Sij(2,2) =            _grad_body_force_z[_qp](2);

  Real SijSij = Sij(0,0) * Sij(0,0) + Sij(0,1) * Sij(1,0) + Sij(0,2) * Sij(2,0) +
                Sij(1,0) * Sij(0,1) + Sij(1,1) * Sij(1,1) + Sij(1,2) * Sij(2,1) +
                Sij(2,0) * Sij(0,2) + Sij(2,1) * Sij(1,2) + Sij(2,2) * Sij(2,2);

  RealTensorValue gij;
  gij(0,0) = _grad_body_force_x[_qp](0);
  gij(0,1) = _grad_body_force_x[_qp](1);
  gij(0,2) = _grad_body_force_x[_qp](2);
  gij(1,0) = _grad_body_force_y[_qp](0);
  gij(1,1) = _grad_body_force_y[_qp](1);
  gij(1,2) = _grad_body_force_y[_qp](2);
  gij(2,0) = _grad_body_force_z[_qp](0);
  gij(2,1) = _grad_body_force_z[_qp](1);
  gij(2,2) = _grad_body_force_z[_qp](2);

  RealTensorValue gij_squared;
  gij_squared(0,0) = gij(0,0) * gij(0,0) + gij(0,1) * gij(1,0) + gij(0,2) * gij(2,0);
  gij_squared(0,1) = gij(0,0) * gij(0,1) + gij(0,1) * gij(1,1) + gij(0,2) * gij(2,1);
  gij_squared(0,2) = gij(0,0) * gij(0,2) + gij(0,1) * gij(1,2) + gij(0,2) * gij(2,2);
  gij_squared(1,0) = gij(1,0) * gij(0,0) + gij(1,1) * gij(1,0) + gij(1,2) * gij(2,0);
  gij_squared(1,1) = gij(1,0) * gij(0,1) + gij(1,1) * gij(1,1) + gij(1,2) * gij(2,1);
  gij_squared(1,2) = gij(1,0) * gij(0,2) + gij(1,1) * gij(1,2) + gij(1,2) * gij(2,2);
  gij_squared(2,0) = gij(2,0) * gij(0,0) + gij(2,1) * gij(1,0) + gij(2,2) * gij(2,0);
  gij_squared(2,1) = gij(2,0) * gij(0,1) + gij(2,1) * gij(1,1) + gij(2,2) * gij(2,1);
  gij_squared(2,2) = gij(2,0) * gij(0,2) + gij(2,1) * gij(1,2) + gij(2,2) * gij(2,2);

  RealTensorValue gji = gij.transpose();

  RealTensorValue gji_squared;
  gji_squared(0,0) = gji(0,0) * gji(0,0) + gji(0,1) * gji(1,0) + gji(0,2) * gji(2,0);
  gji_squared(0,1) = gji(0,0) * gji(0,1) + gji(0,1) * gji(1,1) + gji(0,2) * gji(2,1);
  gji_squared(0,2) = gji(0,0) * gji(0,2) + gji(0,1) * gji(1,2) + gji(0,2) * gji(2,2);
  gji_squared(1,0) = gji(1,0) * gji(0,0) + gji(1,1) * gji(1,0) + gji(1,2) * gji(2,0);
  gji_squared(1,1) = gji(1,0) * gji(0,1) + gji(1,1) * gji(1,1) + gji(1,2) * gji(2,1);
  gji_squared(1,2) = gji(1,0) * gji(0,2) + gji(1,1) * gji(1,2) + gji(1,2) * gji(2,2);
  gji_squared(2,0) = gji(2,0) * gji(0,0) + gji(2,1) * gji(1,0) + gji(2,2) * gji(2,0);
  gji_squared(2,1) = gji(2,0) * gji(0,1) + gji(2,1) * gji(1,1) + gji(2,2) * gji(2,1);
  gji_squared(2,2) = gji(2,0) * gji(0,2) + gji(2,1) * gji(1,2) + gji(2,2) * gji(2,2);

  RealTensorValue isotropic_part;
  isotropic_part(0,0) = isotropic_part(1,1) = isotropic_part(2,2) = (1.0 / 3.0) * (gij(0,0) * gij(0,0) + gij(1,1) * gij(1,1) + gij(2,2) * gij(2,2));
  isotropic_part(0,1) = isotropic_part(1,0) = isotropic_part(0,2) = isotropic_part(2,0) = isotropic_part(1,2) = isotropic_part(2,1) = 0.0;

  RealTensorValue Sdij = 0.5 * (gij_squared + gji_squared) - isotropic_part;

  Real SdijSdij = Sdij(0,0) * Sdij(0,0) + Sdij(0,1) * Sdij(1,0) + Sdij(0,2) * Sdij(2,0) +
                  Sdij(1,0) * Sdij(0,1) + Sdij(1,1) * Sdij(1,1) + Sdij(1,2) * Sdij(2,1) +
                  Sdij(2,0) * Sdij(0,2) + Sdij(2,1) * Sdij(1,2) + Sdij(2,2) * Sdij(2,2);

  Real OP = std::pow(SdijSdij, 1.5) / (std::pow(SijSij, 2.5) + std::pow(SdijSdij, 1.25));

  return _C * std::pow(_d[_qp], 2.0) * OP;
}
