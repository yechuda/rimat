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

#include "TwoPointsMinInverseDistanceDirichletBC.h"

template<>
InputParameters validParams<TwoPointsMinInverseDistanceDirichletBC>()
{
  InputParameters p = validParams<NodalBC>();
  p.addRequiredParam<Real>("first_ref_point_x", "x coordinate of the first refernece point");
  p.addParam<Real>("first_ref_point_y", 0.0, "y coordinate of the first refernece point"); // only required in 2D and 3D
  p.addParam<Real>("first_ref_point_z", 0.0, "z coordinate of the first refernece point"); // only required in 3D
  p.addRequiredParam<Real>("second_ref_point_x", "x coordinate of the second refernece point");
  p.addParam<Real>("second_ref_point_y", 0.0, "y coordinate of the second refernece point"); // only required in 2D and 3D
  p.addParam<Real>("second_ref_point_z", 0.0, "z coordinate of the second refernece point"); // only required in 3D
  p.addRequiredParam<Real>("G0", "inverse offset distance");
  return p;
}


TwoPointsMinInverseDistanceDirichletBC::TwoPointsMinInverseDistanceDirichletBC(const InputParameters & parameters) :
  NodalBC(parameters),
  _first_ref_point_x(getParam<Real>("first_ref_point_x")),
  _first_ref_point_y(getParam<Real>("first_ref_point_y")),
  _first_ref_point_z(getParam<Real>("first_ref_point_z")),
  _second_ref_point_x(getParam<Real>("second_ref_point_x")),
  _second_ref_point_y(getParam<Real>("second_ref_point_y")),
  _second_ref_point_z(getParam<Real>("second_ref_point_z")),
  _G0(getParam<Real>("G0"))
{}

Real
TwoPointsMinInverseDistanceDirichletBC::computeQpResidual()
{
  Real x = (*_current_node)(0);
  Real y = (*_current_node)(1);
  Real z = (*_current_node)(2);

  Real first_geom_d_squared = std::pow(x - _first_ref_point_x, 2.0) + std::pow(y - _first_ref_point_y, 2.0) + std::pow(z - _first_ref_point_z, 2.0);
  Real first_geom_d = std::pow(first_geom_d_squared, 0.5);

  Real second_geom_d_squared = std::pow(x - _second_ref_point_x, 2.0) + std::pow(y - _second_ref_point_y, 2.0) + std::pow(z - _second_ref_point_z, 2.0);
  Real second_geom_d = std::pow(second_geom_d_squared, 0.5);

  Real d = std::min(first_geom_d, second_geom_d);
  Real G = 1.0 / (d + (1.0 / _G0));
  return _u[_qp] - G;
}
