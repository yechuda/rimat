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

#include "FluidLES.h"

template<>
InputParameters validParams<FluidLES>()
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

FluidLES::FluidLES(const InputParameters & parameters) :
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
    _Cs(getParam<Real>("Cs"))
{}

void
FluidLES::computeQpProperties()
{
  _rho[_qp] = _rho_param;

  RealTensorValue Sij;
  Sij(0,0) =            _grad_u[_qp](0);
  Sij(0,1) = Sij(1,0) = 0.5 * (_grad_u[_qp](1) + _grad_v[_qp](0));
  Sij(1,1) =            _grad_v[_qp](1);
  Sij(0,2) = Sij(2,0) = 0.5 * (_grad_u[_qp](2) + _grad_w[_qp](0));
  Sij(1,2) = Sij(2,1) = 0.5 * (_grad_v[_qp](2) + _grad_w[_qp](1));
  Sij(2,2) =            _grad_w[_qp](2);

  Real S_mag_squared = 0.0;

  for (unsigned int n = 0; n < 3; n++)
  {
    for (unsigned int m = 0; m < 3; m++)
    {
      S_mag_squared = S_mag_squared + std::pow(Sij(m,n), 2.0);
    }
  }

  Real S_mag = std::pow(S_mag_squared, 0.5);

  Real h = _current_elem->hmax();

  _mu[_qp] = _mu_mol + _rho_param * std::pow(_Cs, 2.0) * std::pow(h, 2.0) * S_mag;
}
