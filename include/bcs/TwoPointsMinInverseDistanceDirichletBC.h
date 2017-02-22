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

#ifndef TWOPOINTSMININVERSEDISTANCEDIRICHLETBC_H
#define TWOPOINTSMININVERSEDISTANCEDIRICHLETBC_H

#include "NodalBC.h"

class TwoPointsMinInverseDistanceDirichletBC;

template<>
InputParameters validParams<TwoPointsMinInverseDistanceDirichletBC>();

/**
 * Boundary condition of a Dirichlet type
 *
 * Sets the value in the node
 */
class TwoPointsMinInverseDistanceDirichletBC : public NodalBC
{
public:
  TwoPointsMinInverseDistanceDirichletBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  /// The value for this BC
  const Real & _first_ref_point_x;
  const Real & _first_ref_point_y;
  const Real & _first_ref_point_z;
  const Real & _second_ref_point_x;
  const Real & _second_ref_point_y;
  const Real & _second_ref_point_z;
  const Real & _G0;
};

#endif /* TWOPOINTSMININVERSEDISTANCEDIRICHLETBC_H */
