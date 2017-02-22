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

#include "OnePointBoundedInverseDistanceDirichletBC.h"

template<>
InputParameters validParams<OnePointBoundedInverseDistanceDirichletBC>()
{
  InputParameters p = validParams<NodalBC>();
  p.addRequiredParam<Real>("ref_point_x", "x coordinate of the refernece point");
  p.addParam<Real>("ref_point_y", 0.0, "y coordinate of the refernece point"); // only required in 2D and 3D
  p.addParam<Real>("ref_point_z", 0.0, "z coordinate of the refernece point"); // only required in 3D
  p.addRequiredParam<Real>("max_dist", "maximum distance");
  p.addRequiredParam<Real>("G0", "inverse offset distance");
  return p;
}


OnePointBoundedInverseDistanceDirichletBC::OnePointBoundedInverseDistanceDirichletBC(const InputParameters & parameters) :
  NodalBC(parameters),
  _ref_point_x(getParam<Real>("ref_point_x")),
  _ref_point_y(getParam<Real>("ref_point_y")),
  _ref_point_z(getParam<Real>("ref_point_z")),
  _max_dist(getParam<Real>("max_dist")),
  _G0(getParam<Real>("G0"))
{}

Real
OnePointBoundedInverseDistanceDirichletBC::computeQpResidual()
{
  Real x = (*_current_node)(0);
  Real y = (*_current_node)(1);
  Real z = (*_current_node)(2);

  Real geom_d_squared = std::pow(x - _ref_point_x, 2.0) + std::pow(y - _ref_point_y, 2.0) + std::pow(z - _ref_point_z, 2.0);
  Real geom_d = std::pow(geom_d_squared, 0.5);
  Real d = std::min(geom_d, _max_dist);
  Real G = 1.0 / (d + (1.0 / _G0));
  return _u[_qp] - G;
}
