/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "SpalartAllmaras.h"

template<>
InputParameters validParams<SpalartAllmaras>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  params.addRequiredCoupledVar("d", "wall distance");

  // Required parameters
  params.addRequiredParam<Real>("rho", "density");
  params.addRequiredParam<Real>("mu_mol", "molecular dynamic viscosiyty");

  return params;
}



SpalartAllmaras::SpalartAllmaras(const InputParameters & parameters) :
  Kernel(parameters),

  // Coupled old variables
  _u_vel(coupledValueOld("u")),
  _v_vel(coupledValueOld("v")),
  _w_vel(coupledValueOld("w")),
  _d(coupledValueOld("d")),

  // Coupled old gradients
  _grad_u_vel(coupledGradientOld("u")),
  _grad_v_vel(coupledGradientOld("v")),
  _grad_w_vel(coupledGradientOld("w")),

  // Required parameters
  _rho(getParam<Real>("rho")),
  _mu_mol(getParam<Real>("mu_mol")),

  // Old value
  _nu_tilde(valueOld())

{
}



Real SpalartAllmaras::computeQpResidual()
{
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);

  Real first_part = U * _grad_u[_qp] * _test[_i][_qp];

  Real nu_mol = _mu_mol / _rho;
  Real chi = _nu_tilde[_qp] / nu_mol;

  Real ct3 = 1.2;
  Real ct4 = 0.5;

  Real ft2 = ct3 * std::exp(-ct4 * std::pow(chi, 2.0));

  // Real doubleWijWij = (_grad_u_vel[_qp](1) - _grad_v_vel[_qp](0)) * (_grad_v_vel[_qp](0) - _grad_u_vel[_qp](1)) +
  //                     (_grad_u_vel[_qp](2) - _grad_w_vel[_qp](0)) * (_grad_w_vel[_qp](0) - _grad_u_vel[_qp](2)) +
  //                     (_grad_v_vel[_qp](2) - _grad_w_vel[_qp](1)) * (_grad_w_vel[_qp](1) - _grad_v_vel[_qp](2));

  Real doubleWijWij = (_grad_w_vel[_qp](1) - _grad_v_vel[_qp](2)) * (_grad_w_vel[_qp](1) - _grad_v_vel[_qp](2)) +
                      (_grad_u_vel[_qp](2) - _grad_w_vel[_qp](0)) * (_grad_u_vel[_qp](2) - _grad_w_vel[_qp](0)) +
                      (_grad_v_vel[_qp](0) - _grad_u_vel[_qp](1)) * (_grad_v_vel[_qp](0) - _grad_u_vel[_qp](1));

  Real Omega = std::pow(doubleWijWij, 0.5);

  Real cnu1 = 7.1;
  Real fnu1 = std::pow(chi, 3.0) / (std::pow(chi, 3.0) + std::pow(cnu1, 3.0));
  Real fnu2 = 1.0 - (chi / (1.0 + chi * fnu1));
  Real kappa = 0.41;

  Real S_bar = ((_nu_tilde[_qp] * fnu2) / (std::pow(kappa, 2.0) * std::pow(_d[_qp], 2.0)));
  Real c2 = 0.7;
  Real c3 = 0.9;
  Real S_tilde;
  if (S_bar >= -c2 * Omega)
    S_tilde = Omega + S_bar;
  else
    S_tilde = Omega + ((Omega * (Omega * std::pow(c2, 2.0) + c3 * S_bar)) / (Omega * (c3 - 2.0 * c2) - S_bar));

  Real cb1 = 0.1355;
  Real A = cb1 * (1.0 - ft2) * S_tilde;

  Real second_part = -A * _u[_qp] * _test[_i][_qp];

  Real r = std::min(_nu_tilde[_qp] / (S_tilde * std::pow(kappa, 2.0) * std::pow(_d[_qp], 2.0)), 10.0);
  Real cw2 = 0.3;
  Real g = r + cw2 * (std::pow(r, 6.0) - r);
  Real cw3 = 2.0;
  Real cg = (1.0 + std::pow(cw3, 6.0)) / (std::pow(g, 6.0) + std::pow(cw3, 6.0)); // constant for code clarity, not included in original Spalart-Allmaras
  Real fw = g * std::pow(cg, 1.0 / 6.0);
  Real cb2 = 0.622;
  Real sigma = 2.0 / 3.0;
  Real cw1 = cb1 / std::pow(kappa, 2.0) + (1.0 + cb2) / sigma;
  Real B = cw1 * fw - cb1 * ft2 / std::pow(kappa, 2.0);

  Real third_part = B / std::pow(_d[_qp], 2.0) * _u[_qp] * _u[_qp] * _test[_i][_qp];

  Real fourth_part = nu_mol / sigma * _grad_u[_qp] * _grad_test[_i][_qp];

  Real fifth_part = 1.0 / sigma * _u[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

  Real sixth_part = -cb2 / sigma * _grad_u[_qp] * _grad_u[_qp] * _test[_i][_qp];


  return first_part + second_part + third_part + fourth_part + fifth_part + sixth_part;
}



Real SpalartAllmaras::computeQpJacobian()
{
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);

  Real first_part = U * _grad_phi[_j][_qp] * _test[_i][_qp];

  Real nu_mol = _mu_mol / _rho;
  Real chi = _nu_tilde[_qp] / nu_mol;

  Real ct3 = 1.2;
  Real ct4 = 0.5;

  Real ft2 = ct3 * std::exp(-ct4 * std::pow(chi, 2.0));

  // Real doubleWijWij = (_grad_u_vel[_qp](1) - _grad_v_vel[_qp](0)) * (_grad_v_vel[_qp](0) - _grad_u_vel[_qp](1)) +
  //                     (_grad_u_vel[_qp](2) - _grad_w_vel[_qp](0)) * (_grad_w_vel[_qp](0) - _grad_u_vel[_qp](2)) +
  //                     (_grad_v_vel[_qp](2) - _grad_w_vel[_qp](1)) * (_grad_w_vel[_qp](1) - _grad_v_vel[_qp](2));

  Real doubleWijWij = (_grad_w_vel[_qp](1) - _grad_v_vel[_qp](2)) * (_grad_w_vel[_qp](1) - _grad_v_vel[_qp](2)) +
                      (_grad_u_vel[_qp](2) - _grad_w_vel[_qp](0)) * (_grad_u_vel[_qp](2) - _grad_w_vel[_qp](0)) +
                      (_grad_v_vel[_qp](0) - _grad_u_vel[_qp](1)) * (_grad_v_vel[_qp](0) - _grad_u_vel[_qp](1));

  Real Omega = std::pow(doubleWijWij, 0.5);

  Real cnu1 = 7.1;
  Real fnu1 = std::pow(chi, 3.0) / (std::pow(chi, 3.0) + std::pow(cnu1, 3.0));
  Real fnu2 = 1.0 - (chi / (1.0 + chi * fnu1));
  Real kappa = 0.41;

  Real S_bar = ((_nu_tilde[_qp] * fnu2) / (std::pow(kappa, 2.0) * std::pow(_d[_qp], 2.0)));
  Real c2 = 0.7;
  Real c3 = 0.9;
  Real S_tilde;
  if (S_bar >= -c2 * Omega)
    S_tilde = Omega + S_bar;
  else
    S_tilde = Omega + ((Omega * (Omega * std::pow(c2, 2.0) + c3 * S_bar)) / (Omega * (c3 - 2.0 * c2) - S_bar));

  Real cb1 = 0.1355;
  Real A = cb1 * (1.0 - ft2) * S_tilde;

  Real second_part = -A * _phi[_j][_qp] * _test[_i][_qp];

  Real r = std::min(_nu_tilde[_qp] / (S_tilde * std::pow(kappa, 2.0) * std::pow(_d[_qp], 2.0)), 10.0);
  Real cw2 = 0.3;
  Real g = r + cw2 * (std::pow(r, 6.0) - r);
  Real cw3 = 2.0;
  Real cg = (1.0 + std::pow(cw3, 6.0)) / (std::pow(g, 6.0) + std::pow(cw3, 6.0)); // constant for code clarity, not included in original Spalart-Allmaras
  Real fw = g * std::pow(cg, 1.0 / 6.0);
  Real cb2 = 0.622;
  Real sigma = 2.0 / 3.0;
  Real cw1 = cb1 / std::pow(kappa, 2.0) + (1.0 + cb2) / sigma;
  Real B = cw1 * fw - cb1 * ft2 / std::pow(kappa, 2.0);

  Real third_part = B / std::pow(_d[_qp], 2.0) * 2.0 * _u[_qp] * _phi[_j][_qp] * _test[_i][_qp];

  Real fourth_part = nu_mol / sigma * _grad_phi[_j][_qp] * _grad_test[_i][_qp];

  Real fifth_part = 1.0 / sigma * _phi[_j][_qp] * _grad_u[_qp] * _grad_test[_i][_qp] +
                    1.0 / sigma * _u[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];

  Real sixth_part = -cb2 / sigma * 2.0 * _grad_u[_qp] * _grad_phi[_j][_qp] * _test[_i][_qp];


  return first_part + second_part + third_part + fourth_part + fifth_part + sixth_part;
}
