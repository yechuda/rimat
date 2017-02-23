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

#include "SpalartAllmarasAux.h"

template<>
InputParameters validParams<SpalartAllmarasAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables
  params.addRequiredCoupledVar("nu_tilde", "Spalart-Allmaras variable");

  // Required parameters
  params.addRequiredParam<Real>("mu_mol", "molecular dynamic viscosity");
  params.addRequiredParam<Real>("rho", "density");

  return params;
}

SpalartAllmarasAux::SpalartAllmarasAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    // Coupled variables
    _nu_tilde(coupledValue("nu_tilde")),

    // Required parameters
    _mu_mol(getParam<Real>("mu_mol")),
    _rho(getParam<Real>("rho"))

{
}

Real SpalartAllmarasAux::computeValue()
{
  Real nu_mol = _mu_mol / _rho;
  Real chi = _nu_tilde[_qp] / nu_mol;
  Real cnu1 = 7.1;
  Real fnu1 = std::pow(chi, 3.0) / (std::pow(chi, 3.0) + std::pow(cnu1, 3.0));
  return _mu_mol + _rho * _nu_tilde[_qp] * fnu1;
}
