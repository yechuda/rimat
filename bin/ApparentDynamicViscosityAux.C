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

#include "ApparentDynamicViscosityAux.h"

template<>
InputParameters validParams<ApparentDynamicViscosityAux>()
{
  InputParameters params = validParams<AuxKernel>();

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

ApparentDynamicViscosityAux::ApparentDynamicViscosityAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    // Coupled variables
    _grad_u_old(coupledGradientOld("u")),
    _grad_v_old(coupledGradientOld("v")),
    _grad_w_old(coupledGradientOld("w")),

    // Required parameters
    _mu_mol(getParam<Real>("mu_mol")),
    _rho(getParam<Real>("rho")),
    _Cs(getParam<Real>("Cs"))

{
}

Real ApparentDynamicViscosityAux::computeValue()
{
  RealTensorValue Sij;
  Sij(0,0) =            _grad_u_old[qp](0);
  Sij(0,1) = Sij(1,0) = 0.5 * (_grad_u_old[qp](1) + _grad_v_old[qp](0));
  Sij(1,1) =            _grad_v_old[qp](1);
  Sij(0,2) = Sij(2,0) = 0.5 * (_grad_u_old[qp](2) + _grad_w_old[qp](0));
  Sij(1,2) = Sij(2,1) = 0.5 * (_grad_v_old[qp](2) + _grad_w_old[qp](1));
  Sij(2,2) =            _grad_w_old[qp](2);

  return 0.0;
}
