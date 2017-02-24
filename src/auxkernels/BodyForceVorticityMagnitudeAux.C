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

#include "BodyForceVorticityMagnitudeAux.h"

template<>
InputParameters validParams<BodyForceVorticityMagnitudeAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables
  params.addRequiredCoupledVar("body_force_x", "The component of body force in x direction");
  params.addCoupledVar("body_force_y", 0, "The component of body force in y direction"); // only required in 2D and 3D
  params.addCoupledVar("body_force_z", 0, "The component of body force in z direction"); // only required in 3D

  return params;
}

BodyForceVorticityMagnitudeAux::BodyForceVorticityMagnitudeAux(const InputParameters & parameters) :
    AuxKernel(parameters),

    // Coupled variables
    _grad_body_force_x(coupledGradient("body_force_x")),
    _grad_body_force_y(coupledGradient("body_force_y")),
    _grad_body_force_z(coupledGradient("body_force_z"))

{
}

Real BodyForceVorticityMagnitudeAux::computeValue()
{
  Real doubleWijWij = (_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2)) * (_grad_body_force_z[_qp](1) - _grad_body_force_y[_qp](2)) +
                      (_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0)) * (_grad_body_force_x[_qp](2) - _grad_body_force_z[_qp](0)) +
                      (_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1)) * (_grad_body_force_y[_qp](0) - _grad_body_force_x[_qp](1));

  return std::pow(doubleWijWij, 0.5);
}
