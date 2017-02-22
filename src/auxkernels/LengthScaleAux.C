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

#include "LengthScaleAux.h"
#include "Assembly.h"

template<>
InputParameters validParams<LengthScaleAux>()
{
  InputParameters params = validParams<AuxKernel>();

  return params;
}

LengthScaleAux::LengthScaleAux(const InputParameters & parameters) :
    AuxKernel(parameters)

{
}

Real LengthScaleAux::computeValue()
{
  Real vol = _current_elem_volume;
  Real h_current = 2.0 * std::pow(vol, 1.0 / 3.0);

  Real h_neighbor_sum = 0.0;
  unsigned int n_faces = 0;

  for (unsigned int side=0; side<_current_elem->n_sides(); side++)
  {
    if (_current_elem->neighbor(side) != NULL)
      {
        Real neighbor_vol = 0.0;

        n_faces++;

        const Elem * neighbor = _current_elem->neighbor(side);
        unsigned int neighbor_side = neighbor->which_neighbor_am_i(_current_elem);

        _assembly.reinitElemAndNeighbor(_current_elem, side, neighbor, neighbor_side);

        neighbor_vol = _assembly.neighborVolume();

        h_neighbor_sum += 2.0 * std::pow(neighbor_vol, 1.0 / 3.0);
      }
  }

  return (h_current + h_neighbor_sum) / (n_faces + 1.0);
}
