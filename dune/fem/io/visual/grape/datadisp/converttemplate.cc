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
#include <dune/fem/io/file/vtkio.hh>

//==================================================
// Include your typedefs here
//==================================================

// type of discrete function tuple to restore
typedef InTupleType GR_InputType;

template <class GR_GridType,
          class InTupleType>
void process(const GR_GridType& grid,
             const InTupleType& data,
             const double time,
             const int timestep,
             const int myRank,
             const int numProcs)
{
  // do somthing with the data
  // Note, that data is destroyed after.
}

// include main program
#include <dune/fem/io/visual/grape/datadisp/dataconvert.cc>
