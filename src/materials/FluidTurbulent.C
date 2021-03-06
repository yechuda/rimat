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

#include "FluidTurbulent.h"
#include "Assembly.h"

template<>
InputParameters validParams<FluidTurbulent>()
{
  InputParameters params = validParams<Material>();

  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D

  // Required parameters
  params.addRequiredParam<Real>("mu_mol", "molecular dynamic viscosity");
  params.addRequiredParam<Real>("rho", "density");
  params.addRequiredParam<Real>("Cs", "Smagorinsky coefficient");

  return params;
}

FluidTurbulent::FluidTurbulent(const InputParameters & parameters) :
    Material(parameters),

    // Declarations
    _rho(declareProperty<Real>("rho")),
    _mu(declareProperty<Real>("mu")),

    // Gradients
    _grad_u(coupledGradient("u")),
    _grad_v(coupledGradient("v")),
    _grad_w(coupledGradient("w")),

    // Required parameters
    _mu_mol(getParam<Real>("mu_mol")),
    _rho_param(getParam<Real>("rho")),
    _Cs(getParam<Real>("Cs")),
    _current_elem_volume(_assembly.elemVolume())
{}

void
FluidTurbulent::computeQpProperties()
{
  _rho[_qp] = _rho_param;

  RealTensorValue Sij;
  Sij(0,0) =            _grad_u[_qp](0);
  Sij(0,1) = Sij(1,0) = 0.5 * (_grad_u[_qp](1) + _grad_v[_qp](0));
  Sij(1,1) =            _grad_v[_qp](1);
  Sij(0,2) = Sij(2,0) = 0.5 * (_grad_u[_qp](2) + _grad_w[_qp](0));
  Sij(1,2) = Sij(2,1) = 0.5 * (_grad_v[_qp](2) + _grad_w[_qp](1));
  Sij(2,2) =            _grad_w[_qp](2);

  Real SijSij = Sij(0,0) * Sij(0,0) + Sij(0,1) * Sij(1,0) + Sij(0,2) * Sij(2,0) +
                Sij(1,0) * Sij(0,1) + Sij(1,1) * Sij(1,1) + Sij(1,2) * Sij(2,1) +
                Sij(2,0) * Sij(0,2) + Sij(2,1) * Sij(1,2) + Sij(2,2) * Sij(2,2);

  RealTensorValue gij;
  gij(0,0) = _grad_u[_qp](0);
  gij(0,1) = _grad_u[_qp](1);
  gij(0,2) = _grad_u[_qp](2);
  gij(1,0) = _grad_v[_qp](0);
  gij(1,1) = _grad_v[_qp](1);
  gij(1,2) = _grad_v[_qp](2);
  gij(2,0) = _grad_w[_qp](0);
  gij(2,1) = _grad_w[_qp](1);
  gij(2,2) = _grad_w[_qp](2);

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

  Real vol = _current_elem_volume;
  Real h = 2.0 * std::pow(vol, 0.33333333);

  _mu[_qp]  = _mu_mol + _rho_param * 10.6 * std::pow(_Cs, 2.0) * std::pow(h, 2.0) * OP;
}
