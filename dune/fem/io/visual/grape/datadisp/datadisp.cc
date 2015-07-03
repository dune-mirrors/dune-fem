//************************************************************
//
//  (C) written and directed by Robert Kloefkorn
//
//************************************************************

namespace Dune {
/** \addtogroup Visualization
 *
 *  Data written in the GraPE format
 *  using the Dune::DataWriter class
 *  can be either directly visualized during
 *  the computation or can be read separately
 *  with the datadisp utility file. More
 *  detail on how to write data to files
 *  can be found in \ref DiscFuncIO.
 *
 *  If data is to be visualized during a computation
 *  the Dune::DataWriter can also be used by
 *  setting the parameter \b fem.io.grapedisplay
 *  to 1. Each time a file is written to disk,
 *  the data is simultaneously displayed in
 *  GraPE. By pressing the exit button in GraPE
 *  the computation is resumed.
 *
 *  For reading data from disk into GraPE the
 *  programm fem/io/visual/grape/datadisp/datadisp.cc
 *  can be used. An example of its usage can be
 *  found in the file
 *  fem/io/visual/grape/datadisp/programtemplate.cc
 *
 *  In the simplest case the following
 *  program can be used:
 *  \code
 *  typedef IOTupleType GR_InputType;
 *  template <class GrapeDispType,
 *            class GR_GridType,
 *            class DestinationType>
 *  void postProcessing(const GrapeDispType& disp,
 *                      const GR_GridType& grid,
 *                      const double time,
 *                      const DestinationType& Uh)
 *  {}
 *  #include <dune/fem/io/visual/grape/datadisp/datadisp.cc>
 *  \endcode
 *
 *  In the file included in the last line
 *  the function main() is defined.
 *  The resulting executable expects the start
 *  and end index of the data files which are to
 *  be read. The file prefix and the directory
 *  are prescribed through runtime parameters with
 *  the same keys used in the Dune::DataWriter
 *  described in \ref DiscFuncIO.
 *
 *  The type GR_InputType should be a tuple
 *  type holding the types of the discrete functions
 *  to be read.
 *  The \b postProcessing method is called after
 *  each data set is read. The corresponding
 *  discrete functions are passed in the third
 *  parameter; the first parameter holds the
 *  grape handle. Further data, e.g., errors
 *  and derived values, can be added to the
 *  display by using the method
 *  disp.addData(...) as described in the
 *  class Dune::GrapeDataDisplay.
 *  An example which adds the difference between
 *  a discrete solution and an analytical
 *  solution can be found in
 *  fem/io/visual/grape/datadisp/errordisplay.hh.
 *
 **/
}

#include <iostream>
#include <vector>
#include <cassert>
#include <string>

#if HAVE_MPI == 1
#warning "Visualization does not work in parallel"
#endif

#include <dune/common/exceptions.hh>
using namespace Dune;

// include grape visualization
#include <dune/grid/io/visual/grapedatadisplay.hh>
#include <dune/grid/io/visual/combinedgrapedisplay.hh>

// include data reading
#include <dune/fem/io/visual/grape/datadisp/printhelp.cc>
#include <dune/fem/io/visual/grape/datadisp/readiotupledata.cc>
#include <dune/fem/io/visual/grape/datadisp/readioparams.cc>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/function/common/discretefunctionadapter.hh>
#include <dune/fem/misc/mpimanager.hh>

int main(int argc, char **argv)
{
  Fem::MPIManager::initialize(argc,argv);
  try {
    Fem::Parameter::append(argc,argv);
    if (argc < 2)
    {
      print_help(argv[0]);
      return(0);
    }

    if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "-help"))
    {
      print_help(argv[0]);
      return(0);
    }
    return readParameterList(argc,argv);
  }
  catch (Dune::Exception& e)
  {
    std::cerr << e << std::endl;
    return 1;
  }
  return 0;
}
