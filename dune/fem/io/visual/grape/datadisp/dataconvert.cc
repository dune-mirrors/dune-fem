//************************************************************
//
//  (C) written and directed by Robert Kloefkorn
//
//************************************************************
#include <iostream>
#include <vector>
#include <cassert>
#include <string>

#if HAVE_MPI == 1
#warning "Visualization does not work in parallel"
#endif

#include <dune/common/exceptions.hh>
using namespace Dune;

// include definition of grid type
// #include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

// include data reading
#include <dune/fem/io/visual/grape/datadisp/printhelp.cc>

// uses readtuple data instead of readiotupledata.
#include <dune/fem/io/visual/grape/datadisp/readtupledata.cc>
#include <dune/fem/io/visual/grape/datadisp/readioparams.cc>
#include <dune/fem/misc/mpimanager.hh>

int main(int argc, char **argv)
{
  Fem::MPIManager::initialize(argc,argv);
  try {
    Fem::Parameter::append(argc,argv);
#ifdef PARAMETER_APPEND_FUNCTION
    PARAMETER_APPEND_FUNCTION();
#endif

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
    return readParameterList(argc,argv,false);
  }
  catch (Dune::Exception& e)
  {
    std::cerr << e << std::endl;
    return 1;
  }
  return 0;
}
