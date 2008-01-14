//************************************************************
//
//  (C) written and directed by Robert Kloefkorn 
//
//************************************************************
#include <config.h>

///////////////////////////////////////////////////
//
// Include your header defining all necessary types 
//
///////////////////////////////////////////////////
#include "header_name.hh"

// type of discrete function tuple to restore 
typedef IOTupleType GR_InputType;

// preprocessing method fill in some method if Uh 
// should be changed before displaying
template <class GrapeDispType, 
          class GR_GridType,
          class DestinationType>
void postProcessing(const GrapeDispType& disp,
                    const GR_GridType& grid,
                    const double time,
                    const DestinationType& Uh)
{
}

// include main program 
#include <dune/fem/io/visual/grape/datadisp/datadisp.cc>
