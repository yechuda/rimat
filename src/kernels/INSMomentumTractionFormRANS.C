/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "INSMomentumTractionFormRANS.h"

template<>
InputParameters validParams<INSMomentumTractionFormRANS>()
{
  InputParameters params = validParams<INSMomentumBaseRANS>();
  return params;
}

INSMomentumTractionFormRANS::INSMomentumTractionFormRANS(const InputParameters & parameters) :
    INSMomentumBaseRANS(parameters)
{
}

Real
INSMomentumTractionFormRANS::computeQpResidualViscousPart()
{
  // The component'th row (or col, it's symmetric) of the viscous stress tensor
  RealVectorValue tau_row;

  switch (_component)
  {
  case 0:
    tau_row(0) = 2.*_grad_u_vel[_qp](0);                    // 2*du/dx1
    tau_row(1) = _grad_u_vel[_qp](1) + _grad_v_vel[_qp](0); // du/dx2 + dv/dx1
    tau_row(2) = _grad_u_vel[_qp](2) + _grad_w_vel[_qp](0); // du/dx3 + dw/dx1
    break;

  case 1:
    tau_row(0) = _grad_v_vel[_qp](0) + _grad_u_vel[_qp](1); // dv/dx1 + du/dx2
    tau_row(1) = 2.*_grad_v_vel[_qp](1);                    // 2*dv/dx2
    tau_row(2) = _grad_v_vel[_qp](2) + _grad_w_vel[_qp](1); // dv/dx3 + dw/dx2
    break;

  case 2:
    tau_row(0) = _grad_w_vel[_qp](0) + _grad_u_vel[_qp](2); // dw/dx1 + du/dx3
    tau_row(1) = _grad_w_vel[_qp](1) + _grad_v_vel[_qp](2); // dw/dx2 + dv/dx3
    tau_row(2) = 2.*_grad_w_vel[_qp](2);                    // 2*dw/dx3
    break;

  default:
    mooseError("Unrecognized _component requested.");
  }

  // The viscous part, _mu * tau : grad(v)
  return _mu[_qp] * (tau_row * _grad_test[_i][_qp]);
}

Real
INSMomentumTractionFormRANS::computeQpJacobianViscousPart()
{
  // Viscous part, full stress tensor.  The extra contribution comes from the "2"
  // on the diagonal of the viscous stress tensor.
  return _mu[_qp] * (_grad_phi[_j][_qp]             * _grad_test[_i][_qp] +
                _grad_phi[_j][_qp](_component) * _grad_test[_i][_qp](_component));
}

Real
INSMomentumTractionFormRANS::computeQpOffDiagJacobianViscousPart(unsigned jvar)
{
  // In Stokes/Laplacian version, off-diag Jacobian entries wrt u,v,w are zero
  if (jvar == _u_vel_var_number)
    return _mu[_qp] * _grad_phi[_j][_qp](_component) * _grad_test[_i][_qp](0);

  else if (jvar == _v_vel_var_number)
    return _mu[_qp] * _grad_phi[_j][_qp](_component) * _grad_test[_i][_qp](1);

  else if (jvar == _w_vel_var_number)
    return _mu[_qp] * _grad_phi[_j][_qp](_component) * _grad_test[_i][_qp](2);

  else if (jvar == _mu_var_number)
    {
      // The component'th row (or col, it's symmetric) of the viscous stress tensor
      RealVectorValue tau_row;

      switch (_component)
      {
      case 0:
        tau_row(0) = 2.*_grad_u_vel[_qp](0);                    // 2*du/dx1
        tau_row(1) = _grad_u_vel[_qp](1) + _grad_v_vel[_qp](0); // du/dx2 + dv/dx1
        tau_row(2) = _grad_u_vel[_qp](2) + _grad_w_vel[_qp](0); // du/dx3 + dw/dx1
        break;

      case 1:
        tau_row(0) = _grad_v_vel[_qp](0) + _grad_u_vel[_qp](1); // dv/dx1 + du/dx2
        tau_row(1) = 2.*_grad_v_vel[_qp](1);                    // 2*dv/dx2
        tau_row(2) = _grad_v_vel[_qp](2) + _grad_w_vel[_qp](1); // dv/dx3 + dw/dx2
        break;

      case 2:
        tau_row(0) = _grad_w_vel[_qp](0) + _grad_u_vel[_qp](2); // dw/dx1 + du/dx3
        tau_row(1) = _grad_w_vel[_qp](1) + _grad_v_vel[_qp](2); // dw/dx2 + dv/dx3
        tau_row(2) = 2.*_grad_w_vel[_qp](2);                    // 2*dw/dx3
        break;

      default:
        mooseError("Unrecognized _component requested.");
      }

      return _phi[_j][_qp] * (tau_row * _grad_test[_i][_qp]);
    }

  else
    return 0;
}
