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

#include "ElemVolumeAux.h"

template<>
InputParameters validParams<ElemVolumeAux>()
{
  InputParameters params = validParams<AuxKernel>();

  return params;
}

ElemVolumeAux::ElemVolumeAux(const InputParameters & parameters) :
    AuxKernel(parameters)

{
}

Real ElemVolumeAux::computeValue()
{
  return _current_elem_volume;
}
