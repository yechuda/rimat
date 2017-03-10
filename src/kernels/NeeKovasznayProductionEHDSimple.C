/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "NeeKovasznayProductionEHDSimple.h"

template<>
InputParameters validParams<NeeKovasznayProductionEHDSimple>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("body_force_x", "The component of body force in x direction");
  params.addCoupledVar("body_force_y", 0, "The component of body force in y direction"); // only required in 2D and 3D
  params.addCoupledVar("body_force_z", 0, "The component of body force in z direction"); // only required in 3D
  params.addRequiredCoupledVar("d", "wall distance");

  // Required parameters
  params.addRequiredParam<Real>("rho", "density");
  params.addRequiredParam<Real>("C", "EHD production coefficient");

  return params;
}



NeeKovasznayProductionEHDSimple::NeeKovasznayProductionEHDSimple(const InputParameters & parameters) :
  Kernel(parameters),

  // Coupled variables
  _d(coupledValue("d")),

  // Coupled gradients
  _grad_body_force_x(coupledGradient("body_force_x")),
  _grad_body_force_y(coupledGradient("body_force_y")),
  _grad_body_force_z(coupledGradient("body_force_z")),

  // Variable numberings
  _body_force_x_var(coupled("body_force_x")),
  _body_force_y_var(coupled("body_force_y")),
  _body_force_z_var(coupled("body_force_z")),
  _d_var(coupled("d")),

  // Required parameters
  _rho(getParam<Real>("rho")),
  _C(getParam<Real>("C"))

{
}



Real NeeKovasznayProductionEHDSimple::computeQpResidual()
{
  Real V_squared = std::pow(_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2), 2.0) + std::pow(_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0), 2.0) + std::pow(_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1), 2.0);
  Real V = std::pow(V_squared, 0.5);

  // return -_C / _rho * std::pow(_d[_qp], 2.0) * V * _test[_i][_qp];
  // return _C / _rho * std::pow(_d[_qp], 2.0) * V * _test[_i][_qp];
  return -_C * V * _test[_i][_qp];
}



Real NeeKovasznayProductionEHDSimple::computeQpJacobian()
{
  return 0.0;
}

Real NeeKovasznayProductionEHDSimple::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _body_force_x_var)
  {
    Real V_squared = std::pow(_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2), 2.0) + std::pow(_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0), 2.0) + std::pow(_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1), 2.0);
    Real V = std::pow(V_squared, 0.5);
    RealVectorValue grad_B_column(_grad_body_force_x[_qp](0), _grad_body_force_y[_qp](0), _grad_body_force_z[_qp](0));
    Real dV_dBxj = (_grad_body_force_x[_qp] - grad_B_column) / V * _grad_phi[_j][_qp];

    // return -_C / _rho * std::pow(_d[_qp], 2.0) * dV_dBxj * _test[_i][_qp];
    // return _C / _rho * std::pow(_d[_qp], 2.0) * dV_dBxj * _test[_i][_qp];
    return -_C * dV_dBxj * _test[_i][_qp];
  }

  else if (jvar == _body_force_y_var)
  {
    Real V_squared = std::pow(_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2), 2.0) + std::pow(_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0), 2.0) + std::pow(_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1), 2.0);
    Real V = std::pow(V_squared, 0.5);
    RealVectorValue grad_B_column(_grad_body_force_x[_qp](1), _grad_body_force_y[_qp](1), _grad_body_force_z[_qp](1));
    Real dV_dByj = (_grad_body_force_y[_qp] - grad_B_column) / V * _grad_phi[_j][_qp];

    // return -_C / _rho * std::pow(_d[_qp], 2.0) * dV_dByj * _test[_i][_qp];
    // return _C / _rho * std::pow(_d[_qp], 2.0) * dV_dByj * _test[_i][_qp];
    return -_C * dV_dByj * _test[_i][_qp];
  }

  else if (jvar == _body_force_z_var)
  {
    Real V_squared = std::pow(_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2), 2.0) + std::pow(_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0), 2.0) + std::pow(_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1), 2.0);
    Real V = std::pow(V_squared, 0.5);
    RealVectorValue grad_B_column(_grad_body_force_x[_qp](2), _grad_body_force_y[_qp](2), _grad_body_force_z[_qp](2));
    Real dV_dBzj = (_grad_body_force_z[_qp] - grad_B_column) / V * _grad_phi[_j][_qp];

    // return -_C / _rho * std::pow(_d[_qp], 2.0) * dV_dBzj * _test[_i][_qp];
    // return _C / _rho * std::pow(_d[_qp], 2.0) * dV_dBzj * _test[_i][_qp];
    return -_C * dV_dBzj * _test[_i][_qp];
  }

  else if (jvar == _d_var)
  {
    // Real V_squared = std::pow(_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2), 2.0) + std::pow(_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0), 2.0) + std::pow(_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1), 2.0);
    // Real V = std::pow(V_squared, 0.5);

    // return -2.0 * _C / _rho * _d[_qp] * _phi[_j][_qp] * V * _test[_i][_qp];
    // return 2.0 * _C / _rho * _d[_qp] * _phi[_j][_qp] * V * _test[_i][_qp];
    return 0.0;
  }

  else
    return 0.0;
}
