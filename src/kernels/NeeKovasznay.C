/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "NeeKovasznay.h"

template<>
InputParameters validParams<NeeKovasznay>()
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
  params.addRequiredParam<Real>("A", "production coefficient");
  params.addRequiredParam<Real>("B", "destruction coefficient");

  return params;
}



NeeKovasznay::NeeKovasznay(const InputParameters & parameters) :
  Kernel(parameters),

  // Coupled variables
  _u_vel(coupledValue("u")),
  _v_vel(coupledValue("v")),
  _w_vel(coupledValue("w")),
  _d(coupledValue("d")),

  // Coupled gradients
  _grad_u_vel(coupledGradient("u")),
  _grad_v_vel(coupledGradient("v")),
  _grad_w_vel(coupledGradient("w")),

  // Variable numberings
  _u_vel_var_number(coupled("u")),
  _v_vel_var_number(coupled("v")),
  _w_vel_var_number(coupled("w")),
  _d_var_number(coupled("d")),

  // Required parameters
  _rho(getParam<Real>("rho")),
  _mu_mol(getParam<Real>("mu_mol")),
  _A(getParam<Real>("A")),
  _B(getParam<Real>("B"))

{
}



Real NeeKovasznay::computeQpResidual()
{
  // convection part
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  Real convection_part = U * _grad_u[_qp] * _test[_i][_qp];

  // diffusion part
  Real diffusion_part = _u[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

  // production part
  Real S_squared = std::pow(_grad_w_vel[_qp](1) - _grad_v_vel[_qp](2), 2.0) + std::pow(_grad_u_vel[_qp](2) - _grad_w_vel[_qp](0), 2.0) + std::pow(_grad_v_vel[_qp](0) - _grad_u_vel[_qp](1), 2.0);
  Real S = std::pow(S_squared, 0.5);
  Real nu_mol = _mu_mol / _rho;
  Real production_part = -_A * (_u[_qp] - nu_mol) * S * _test[_i][_qp];

  // destruction part
  Real destruction_part = _B / std::pow(_d[_qp], 2.0) * _u[_qp] * (_u[_qp] - nu_mol) * _test[_i][_qp];

  return convection_part + diffusion_part + production_part + destruction_part;
}



Real NeeKovasznay::computeQpJacobian()
{
  // convection part
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  Real convection_part = U * _grad_phi[_j][_qp] * _test[_i][_qp];

  // diffusion part
  Real diffusion_part = _phi[_j][_qp] * _grad_u[_qp] * _grad_test[_i][_qp] +
                        _u[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp];

  // production part
  Real S_squared = std::pow(_grad_w_vel[_qp](1) - _grad_v_vel[_qp](2), 2.0) + std::pow(_grad_u_vel[_qp](2) - _grad_w_vel[_qp](0), 2.0) + std::pow(_grad_v_vel[_qp](0) - _grad_u_vel[_qp](1), 2.0);
  Real S = std::pow(S_squared, 0.5);
  Real nu_mol = _mu_mol / _rho;
  Real production_part = -_A * _phi[_j][_qp] * S * _test[_i][_qp];

  // destruction part
  Real destruction_part = _B / std::pow(_d[_qp], 2.0) * 2.0 * _u[_qp] * _phi[_j][_qp] * _test[_i][_qp] -
                          _B / std::pow(_d[_qp], 2.0) * _phi[_j][_qp] * nu_mol * _test[_i][_qp];

  return convection_part + diffusion_part + production_part + destruction_part;
}

Real NeeKovasznay::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _u_vel_var_number)
  {
    // convection part
    Real convection_part = _phi[_j][_qp] * _grad_u[_qp](0) * _test[_i][_qp];

    // production part
    Real S_squared = std::pow(_grad_w_vel[_qp](1) - _grad_v_vel[_qp](2), 2.0) + std::pow(_grad_u_vel[_qp](2) - _grad_w_vel[_qp](0), 2.0) + std::pow(_grad_v_vel[_qp](0) - _grad_u_vel[_qp](1), 2.0);
    Real S = std::pow(S_squared, 0.5);
    Real nu_mol = _mu_mol / _rho;
    RealVectorValue grad_U_column(_grad_u_vel[_qp](0), _grad_v_vel[_qp](0), _grad_w_vel[_qp](0));
    Real dS_duj = ((_grad_u_vel[_qp] - grad_U_column) / S) * _grad_phi[_j][_qp];
    Real production_part = -_A * (_u[_qp] - nu_mol) * dS_duj * _test[_i][_qp];

    return convection_part + production_part;
  }

  else if (jvar == _v_vel_var_number)
  {
    // convection part
    Real convection_part = _phi[_j][_qp] * _grad_u[_qp](1) * _test[_i][_qp];

    // production part
    Real S_squared = std::pow(_grad_w_vel[_qp](1) - _grad_v_vel[_qp](2), 2.0) + std::pow(_grad_u_vel[_qp](2) - _grad_w_vel[_qp](0), 2.0) + std::pow(_grad_v_vel[_qp](0) - _grad_u_vel[_qp](1), 2.0);
    Real S = std::pow(S_squared, 0.5);
    Real nu_mol = _mu_mol / _rho;
    RealVectorValue grad_U_column(_grad_u_vel[_qp](1), _grad_v_vel[_qp](1), _grad_w_vel[_qp](1));
    Real dS_dvj = ((_grad_v_vel[_qp] - grad_U_column) / S) * _grad_phi[_j][_qp];
    Real production_part = -_A * (_u[_qp] - nu_mol) * dS_dvj * _test[_i][_qp];

    return convection_part + production_part;
  }

  else if (jvar == _w_vel_var_number)
  {
    // convection part
    Real convection_part = _phi[_j][_qp] * _grad_u[_qp](2) * _test[_i][_qp];

    // production part
    Real S_squared = std::pow(_grad_w_vel[_qp](1) - _grad_v_vel[_qp](2), 2.0) + std::pow(_grad_u_vel[_qp](2) - _grad_w_vel[_qp](0), 2.0) + std::pow(_grad_v_vel[_qp](0) - _grad_u_vel[_qp](1), 2.0);
    Real S = std::pow(S_squared, 0.5);
    Real nu_mol = _mu_mol / _rho;
    RealVectorValue grad_U_column(_grad_u_vel[_qp](2), _grad_v_vel[_qp](2), _grad_w_vel[_qp](2));
    Real dS_dwj = ((_grad_w_vel[_qp] - grad_U_column) / S) * _grad_phi[_j][_qp];
    Real production_part = -_A * (_u[_qp] - nu_mol) * dS_dwj * _test[_i][_qp];

    return convection_part + production_part;
  }

  else if (jvar == _d_var_number)
  {
    // destruction part
    Real nu_mol = _mu_mol / _rho;
    Real destruction_part = -2.0 * _B / std::pow(_d[_qp], 3.0) * _phi[_j][_qp] * _u[_qp] * (_u[_qp] - nu_mol) * _test[_i][_qp];

    return destruction_part;
  }

  else
    return 0.0;
}
