/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "INSkTransport.h"

template<>
InputParameters validParams<INSkTransport>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  params.addRequiredCoupledVar("epsilon", "dissipation rate");

  // Required parameters
  params.addParam<Real>("Cmu", 0.09, "Cmu coefficient");

  return params;
}



INSkTransport::INSkTransport(const InputParameters & parameters) :
  Kernel(parameters),

  // Coupled variables
  _u_vel(coupledValue("u")),
  _v_vel(coupledValue("v")),
  _w_vel(coupledValue("w")),
  _epsilon(coupledValue("epsilon")),

  // Gradients
  _grad_u_vel(coupledGradient("u")),
  _grad_v_vel(coupledGradient("v")),
  _grad_w_vel(coupledGradient("w")),

  // Variable numberings
  _u_vel_var_number(coupled("u")),
  _v_vel_var_number(coupled("v")),
  _w_vel_var_number(coupled("w")),
  _epsilon_var_number(coupled("epsilon")),

  // Required parameters
  _Cmu(getParam<Real>("Cmu"))

{
}



Real INSkTransport::computeQpResidual()
{
  // convection part
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  Real convection_part = U * _grad_u[_qp] * _test[_i][_qp];

  // diffusion part
  Real diffusion_part = _Cmu / _epsilon[_qp] * _u[_qp] * _u[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

  // production part
  Real S = 2.0 * _grad_u_vel[_qp](0) * _grad_u_vel[_qp](0) + _grad_u_vel[_qp](1) * (_grad_v_vel[_qp](0) + _grad_u_vel[_qp](1)) + _grad_u_vel[_qp](2) * (_grad_w_vel[_qp](0) + _grad_u_vel[_qp](2)) +
           _grad_v_vel[_qp](0) * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0)) + 2.0 * _grad_v_vel[_qp](1) * _grad_v_vel[_qp](1) + _grad_v_vel[_qp](2) * (_grad_w_vel[_qp](1) + _grad_v_vel[_qp](2)) +
           _grad_w_vel[_qp](0) * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0)) + _grad_w_vel[_qp](1) * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1)) + 2.0 * _grad_w_vel[_qp](2) * _grad_w_vel[_qp](2);

  Real production_part = -_Cmu / _epsilon[_qp] * _u[_qp] * _u[_qp] * S * _test[_i][_qp];

  // destruction part
  Real destruction_part = _epsilon[_qp] * _test[_i][_qp];

  return convection_part + diffusion_part + production_part + destruction_part;
}


Real INSkTransport::computeQpJacobian()
{
  // convection part
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  Real convection_part = U * _grad_phi[_j][_qp] * _test[_i][_qp];

  // diffusion part
  Real diffusion_part = _Cmu / _epsilon[_qp] * _u[_qp] * _u[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp] +
                        _Cmu / _epsilon[_qp] * 2.0 * _u[_qp] * _phi[_j][_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

  // production part
  Real S = 2.0 * _grad_u_vel[_qp](0) * _grad_u_vel[_qp](0) + _grad_u_vel[_qp](1) * (_grad_v_vel[_qp](0) + _grad_u_vel[_qp](1)) + _grad_u_vel[_qp](2) * (_grad_w_vel[_qp](0) + _grad_u_vel[_qp](2)) +
           _grad_v_vel[_qp](0) * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0)) + 2.0 * _grad_v_vel[_qp](1) * _grad_v_vel[_qp](1) + _grad_v_vel[_qp](2) * (_grad_w_vel[_qp](1) + _grad_v_vel[_qp](2)) +
           _grad_w_vel[_qp](0) * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0)) + _grad_w_vel[_qp](1) * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1)) + 2.0 * _grad_w_vel[_qp](2) * _grad_w_vel[_qp](2);

  Real production_part = -_Cmu / _epsilon[_qp] * 2.0 * _u[_qp] * _phi[_j][_qp] * S * _test[_i][_qp];

  return convection_part + diffusion_part + production_part;
}


Real INSkTransport::computeQpOffDiagJacobian(unsigned jvar)
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

    Real production_part = -_Cmu / _epsilon[_qp] * _u[_qp] * _u[_qp] * 2.0 * _grad_phi[_j][_qp] * tau_row * _test[_i][_qp];

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

    Real production_part = -_Cmu / _epsilon[_qp] * _u[_qp] * _u[_qp] * 2.0 * _grad_phi[_j][_qp] * tau_row * _test[_i][_qp];

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

    Real production_part = -_Cmu / _epsilon[_qp] * _u[_qp] * _u[_qp] * 2.0 * _grad_phi[_j][_qp] * tau_row * _test[_i][_qp];

    return convection_part + production_part;
  }

  else if (jvar == _epsilon_var_number)
  {
    // diffusion part
    Real diffusion_part = -_Cmu / (_epsilon[_qp] * _epsilon[_qp]) * _phi[_j][_qp] * _u[_qp] * _u[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

    // production part
    Real S = 2.0 * _grad_u_vel[_qp](0) * _grad_u_vel[_qp](0) + _grad_u_vel[_qp](1) * (_grad_v_vel[_qp](0) + _grad_u_vel[_qp](1)) + _grad_u_vel[_qp](2) * (_grad_w_vel[_qp](0) + _grad_u_vel[_qp](2)) +
             _grad_v_vel[_qp](0) * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0)) + 2.0 * _grad_v_vel[_qp](1) * _grad_v_vel[_qp](1) + _grad_v_vel[_qp](2) * (_grad_w_vel[_qp](1) + _grad_v_vel[_qp](2)) +
             _grad_w_vel[_qp](0) * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0)) + _grad_w_vel[_qp](1) * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1)) + 2.0 * _grad_w_vel[_qp](2) * _grad_w_vel[_qp](2);

    Real production_part = _Cmu / (_epsilon[_qp] * _epsilon[_qp]) * _phi[_j][_qp] * _u[_qp] * _u[_qp] * S * _test[_i][_qp];

    // destruction part
    Real destruction_part = _phi[_j][_qp] * _test[_i][_qp];

    return diffusion_part + production_part + destruction_part;
  }

  else
    return 0.0;
}
