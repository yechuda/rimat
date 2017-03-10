/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "NeeKovasznayProductionEHDWALE.h"

template<>
InputParameters validParams<NeeKovasznayProductionEHDWALE>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled variables
  params.addRequiredCoupledVar("body_force_x", "The component of body force in x direction");
  params.addCoupledVar("body_force_y", 0, "The component of body force in y direction"); // only required in 2D and 3D
  params.addCoupledVar("body_force_z", 0, "The component of body force in z direction"); // only required in 3D
  params.addRequiredCoupledVar("d", "wall distance");

  // Required parameters
  params.addRequiredParam<Real>("C", "EHD production coefficient");

  return params;
}



NeeKovasznayProductionEHDWALE::NeeKovasznayProductionEHDWALE(const InputParameters & parameters) :
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
  _C(getParam<Real>("C"))

{
}



Real NeeKovasznayProductionEHDWALE::computeQpResidual()
{
  RealTensorValue Sij;
  Sij(0,0) =            _grad_body_force_x[_qp](0);
  Sij(0,1) = Sij(1,0) = 0.5 * (_grad_body_force_x[_qp](1) + _grad_body_force_y[_qp](0));
  Sij(1,1) =            _grad_body_force_y[_qp](1);
  Sij(0,2) = Sij(2,0) = 0.5 * (_grad_body_force_x[_qp](2) + _grad_body_force_z[_qp](0));
  Sij(1,2) = Sij(2,1) = 0.5 * (_grad_body_force_y[_qp](2) + _grad_body_force_z[_qp](1));
  Sij(2,2) =            _grad_body_force_z[_qp](2);

  Real SijSij = Sij(0,0) * Sij(0,0) + Sij(0,1) * Sij(1,0) + Sij(0,2) * Sij(2,0) +
                Sij(1,0) * Sij(0,1) + Sij(1,1) * Sij(1,1) + Sij(1,2) * Sij(2,1) +
                Sij(2,0) * Sij(0,2) + Sij(2,1) * Sij(1,2) + Sij(2,2) * Sij(2,2);

  RealTensorValue gij;
  gij(0,0) = _grad_body_force_x[_qp](0);
  gij(0,1) = _grad_body_force_x[_qp](1);
  gij(0,2) = _grad_body_force_x[_qp](2);
  gij(1,0) = _grad_body_force_y[_qp](0);
  gij(1,1) = _grad_body_force_y[_qp](1);
  gij(1,2) = _grad_body_force_y[_qp](2);
  gij(2,0) = _grad_body_force_z[_qp](0);
  gij(2,1) = _grad_body_force_z[_qp](1);
  gij(2,2) = _grad_body_force_z[_qp](2);

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

  return -_C * std::pow(_d[_qp], 2.0) * OP * _test[_i][_qp];
}



Real NeeKovasznayProductionEHDWALE::computeQpJacobian()
{
  return 0.0;
}

Real NeeKovasznayProductionEHDWALE::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _body_force_x_var)
  {
    // Real V_squared = std::pow(_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2), 2.0) + std::pow(_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0), 2.0) + std::pow(_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1), 2.0);
    // Real V = std::pow(V_squared, 0.5);
    // RealVectorValue grad_B_column(_grad_body_force_x[_qp](0), _grad_body_force_y[_qp](0), _grad_body_force_z[_qp](0));
    // Real dV_dBxj = (_grad_body_force_x[_qp] - grad_B_column) / V * _grad_phi[_j][_qp];
    //
    // return -_C / _rho * std::pow(_d[_qp], 2.0) * dV_dBxj * _test[_i][_qp];

    return 0.0;
  }

  else if (jvar == _body_force_y_var)
  {
    // Real V_squared = std::pow(_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2), 2.0) + std::pow(_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0), 2.0) + std::pow(_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1), 2.0);
    // Real V = std::pow(V_squared, 0.5);
    // RealVectorValue grad_B_column(_grad_body_force_x[_qp](1), _grad_body_force_y[_qp](1), _grad_body_force_z[_qp](1));
    // Real dV_dByj = (_grad_body_force_y[_qp] - grad_B_column) / V * _grad_phi[_j][_qp];
    //
    // return -_C / _rho * std::pow(_d[_qp], 2.0) * dV_dByj * _test[_i][_qp];

    return 0.0;
  }

  else if (jvar == _body_force_z_var)
  {
    // Real V_squared = std::pow(_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2), 2.0) + std::pow(_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0), 2.0) + std::pow(_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1), 2.0);
    // Real V = std::pow(V_squared, 0.5);
    // RealVectorValue grad_B_column(_grad_body_force_x[_qp](2), _grad_body_force_y[_qp](2), _grad_body_force_z[_qp](2));
    // Real dV_dBzj = (_grad_body_force_z[_qp] - grad_B_column) / V * _grad_phi[_j][_qp];
    //
    // return -_C / _rho * std::pow(_d[_qp], 2.0) * dV_dBzj * _test[_i][_qp];

    return 0.0;
  }

  else if (jvar == _d_var)
  {
    // Real V_squared = std::pow(_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2), 2.0) + std::pow(_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0), 2.0) + std::pow(_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1), 2.0);
    // Real V = std::pow(V_squared, 0.5);
    //
    // return -2.0 * _C / _rho * _d[_qp] * _phi[_j][_qp] * V * _test[_i][_qp];

    return 0.0;
  }

  else
    return 0.0;
}
