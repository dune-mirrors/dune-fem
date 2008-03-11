//************************************************************
  // Example
  // ./datadisp 0 5 -m solution -p . paramfile:parameter 
//************************************************************
#include <config.h>


///////////////////////////////////////////////////
//
// Include your header defining all necessary types 
//
///////////////////////////////////////////////////
#include "models.hh"
// type of discrete function tuple to restore 
typedef IOTupleType GR_InputType;

#include <dune/fem/io/visual/grape/datadisp/errordisplay.hh>
// preprocessing method fill in some method if Uh 
// should be changed before displaying
template <class GrapeDispType, 
          class GR_GridType,
          class DestinationType>
void postProcessing(GrapeDispType& disp,
                    const GR_GridType& grid,
                    const double time,
                    const DestinationType& Uh)
{
  DisplayErrorFunction<DestinationType,InitialDataType>::
    apply(disp,grid,time,Uh);
}
// include main program 
#include <dune/fem/io/visual/grape/datadisp/datadisp.cc>

