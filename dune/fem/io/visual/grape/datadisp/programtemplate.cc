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
// An example implementation computing the error of an
// approximation can be found in the apply method of the class
// DisplayErrorFunction in file errordisplay.hh
template <class GrapeDispType,
          class GR_GridType>
void postProcessing(const GrapeDispType& disp,
                    const GR_GridType& grid,
                    const double time,
                    const IOTupleType& data)
{
}

// include main program
#include <dune/fem/io/visual/grape/datadisp/datadisp.cc>
