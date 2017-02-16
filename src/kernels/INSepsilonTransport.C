/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "INSepsilonTransport.h"

template<>
InputParameters validParams<INSepsilonTransport>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  params.addRequiredCoupledVar("k", "mean turbulent kinetic energy");

  // Required parameters
  params.addParam<Real>("Cmu", 0.09, "Cmu coefficient");
  params.addParam<Real>("sigma_epsilon", 1.3, "sigma_epsilon coefficient");
  params.addParam<Real>("C_epsilon1", 1.44, "C_epsilon1 coefficient");
  params.addParam<Real>("C_epsilon2", 1.92, "C_epsilon2 coefficient");

  return params;
}



INSepsilonTransport::INSepsilonTransport(const InputParameters & parameters) :
  Kernel(parameters),

  // Coupled variables
  _u_vel(coupledValue("u")),
  _v_vel(coupledValue("v")),
  _w_vel(coupledValue("w")),
  _k(coupledValue("k")),

  // Gradients
  _grad_u_vel(coupledGradient("u")),
  _grad_v_vel(coupledGradient("v")),
  _grad_w_vel(coupledGradient("w")),

  // Variable numberings
  _u_vel_var_number(coupled("u")),
  _v_vel_var_number(coupled("v")),
  _w_vel_var_number(coupled("w")),
  _k_var_number(coupled("k")),

  // Required parameters
  _Cmu(getParam<Real>("Cmu")),
  _sigma_epsilon(getParam<Real>("sigma_epsilon")),
  _C_epsilon1(getParam<Real>("C_epsilon1")),
  _C_epsilon2(getParam<Real>("C_epsilon2"))

{
}



Real INSepsilonTransport::computeQpResidual()
{
  // convection part
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  Real convection_part = U * _grad_u[_qp] * _test[_i][_qp];

  // diffusion part
  Real diffusion_part = (_Cmu / _sigma_epsilon) * _k[_qp] * _k[_qp] / _u[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

  // production part
  Real S = 2.0 * _grad_u_vel[_qp](0) * _grad_u_vel[_qp](0) + _grad_u_vel[_qp](1) * (_grad_v_vel[_qp](0) + _grad_u_vel[_qp](1)) + _grad_u_vel[_qp](2) * (_grad_w_vel[_qp](0) + _grad_u_vel[_qp](2)) +
          _grad_v_vel[_qp](0) * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0)) + 2.0 * _grad_v_vel[_qp](1) * _grad_v_vel[_qp](1) + _grad_v_vel[_qp](2) * (_grad_w_vel[_qp](1) + _grad_v_vel[_qp](2)) +
          _grad_w_vel[_qp](0) * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0)) + _grad_w_vel[_qp](1) * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1)) + 2.0 * _grad_w_vel[_qp](2) * _grad_w_vel[_qp](2);

  Real production_part = -_Cmu * _C_epsilon1 * _k[_qp] * S * _test[_i][_qp];

  // destruction part
  Real destruction_part = _C_epsilon2 / _k[_qp] * _u[_qp] * _u[_qp] * _test[_i][_qp];

  return convection_part + diffusion_part + production_part + destruction_part;
}


Real INSepsilonTransport::computeQpJacobian()
{
  // convection part
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  Real convection_part = U * _grad_phi[_j][_qp] * _test[_i][_qp];

  // diffusion part
  Real diffusion_part = (_Cmu / _sigma_epsilon) * _k[_qp] * _k[_qp] / _u[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp] - (_Cmu / _sigma_epsilon) * _k[_qp] * _k[_qp] / (_u[_qp] * _u[_qp]) * _phi[_j][_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

  // destruction part
  Real destruction_part = _C_epsilon2 / _k[_qp] * 2.0 * _u[_qp] * _phi[_j][_qp] * _test[_i][_qp];

  return convection_part + diffusion_part + destruction_part;
}


Real INSepsilonTransport::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _u_vel_var_number)
  {
    // convection part
    Real convection_part = _phi[_j][_qp] * _grad_u[_qp](0) * _test[_i][_qp];

    // production part
    RealVectorValue tau_row;
    tau_row(0) = 2.0 * _grad_u_vel[_qp](0);                 // 2*du/dx1
    tau_row(1) = _grad_u_vel[_qp](1) + _grad_v_vel[_qp](0); // du/dx2 + dv/dx1
    tau_row(2) = _grad_u_vel[_qp](2) + _grad_w_vel[_qp](0); // du/dx3 + dw/dx1

    Real production_part = -_Cmu * _C_epsilon1 * _k[_qp] * 2.0 * _grad_phi[_j][_qp] * tau_row * _test[_i][_qp];

    return convection_part + production_part;
  }

  else if (jvar == _v_vel_var_number)
  {
    // convection part
    Real convection_part = _phi[_j][_qp] * _grad_u[_qp](1) * _test[_i][_qp];

    // production part
    RealVectorValue tau_row;
    tau_row(0) = _grad_v_vel[_qp](0) + _grad_u_vel[_qp](1); // dv/dx1 + du/dx2
    tau_row(1) = 2.0 * _grad_v_vel[_qp](1);                 // 2*dv/dx2
    tau_row(2) = _grad_v_vel[_qp](2) + _grad_w_vel[_qp](1); // dv/dx3 + dw/dx2

    Real production_part = -_Cmu * _C_epsilon1 * _k[_qp] * 2.0 * _grad_phi[_j][_qp] * tau_row * _test[_i][_qp];

    return convection_part + production_part;
  }

  else if (jvar == _w_vel_var_number)
  {
    // convection part
    Real convection_part = _phi[_j][_qp] * _grad_u[_qp](2) * _test[_i][_qp];

    // production part
    RealVectorValue tau_row;
    tau_row(0) = _grad_w_vel[_qp](0) + _grad_u_vel[_qp](2); // dw/dx1 + du/dx3
    tau_row(1) = _grad_w_vel[_qp](1) + _grad_v_vel[_qp](2); // dw/dx2 + dv/dx3
    tau_row(2) = 2.0 * _grad_w_vel[_qp](2);                 // 2*dw/dx3

    Real production_part = -_Cmu * _C_epsilon1 * _k[_qp] * 2.0 * _grad_phi[_j][_qp] * tau_row * _test[_i][_qp];

    return convection_part + production_part;
  }

  else if (jvar == _k_var_number)
  {
    // diffusion part
    Real diffusion_part = (_Cmu / _sigma_epsilon) * 2.0 * _k[_qp] * _phi[_j][_qp] / _u[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

    // production part
    Real S = 2.0 * _grad_u_vel[_qp](0) * _grad_u_vel[_qp](0) + _grad_u_vel[_qp](1) * (_grad_v_vel[_qp](0) + _grad_u_vel[_qp](1)) + _grad_u_vel[_qp](2) * (_grad_w_vel[_qp](0) + _grad_u_vel[_qp](2)) +
            _grad_v_vel[_qp](0) * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0)) + 2.0 * _grad_v_vel[_qp](1) * _grad_v_vel[_qp](1) + _grad_v_vel[_qp](2) * (_grad_w_vel[_qp](1) + _grad_v_vel[_qp](2)) +
            _grad_w_vel[_qp](0) * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0)) + _grad_w_vel[_qp](1) * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1)) + 2.0 * _grad_w_vel[_qp](2) * _grad_w_vel[_qp](2);

    Real production_part = -_Cmu * _C_epsilon1 * _phi[_j][_qp] * S * _test[_i][_qp];

    // destruction part
    Real destruction_part = -_C_epsilon2 / (_k[_qp] * _k[_qp]) * _phi[_j][_qp] * _u[_qp] * _u[_qp] * _test[_i][_qp];

    return diffusion_part + production_part + destruction_part;
  }

  else
    return 0.0;
}
