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

#ifndef ONEPOINTBOUNDEDINVERSEDISTANCEDIRICHLETBC_H
#define ONEPOINTBOUNDEDINVERSEDISTANCEDIRICHLETBC_H

#include "NodalBC.h"

class OnePointBoundedInverseDistanceDirichletBC;

template<>
InputParameters validParams<OnePointBoundedInverseDistanceDirichletBC>();

/**
 * Boundary condition of a Dirichlet type
 *
 * Sets the value in the node
 */
class OnePointBoundedInverseDistanceDirichletBC : public NodalBC
{
public:
  OnePointBoundedInverseDistanceDirichletBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  /// The value for this BC
  const Real & _ref_point_x;
  const Real & _ref_point_y;
  const Real & _ref_point_z;
  const Real & _max_dist;
  const Real & _G0;
};

#endif /* ONEPOINTBOUNDEDINVERSEDISTANCEDIRICHLETBC_H */
