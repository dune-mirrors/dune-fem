// undef NDEBUG so we can always use assert.
#undef NDEBUG

#ifdef YASPGRID
#define LEAFGRIDPART
#endif

#include <config.h>
#include <iostream>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>

#include <dune/fem/misc/double.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/io/streams/binarystreams.hh>
#include <dune/fem/io/streams/asciistreams.hh>
#include <dune/fem/io/streams/virtualstreams.hh>
#include <dune/fem/io/file/vtkio.hh>

#include "testgrid.hh"
#include "exactsolution.hh"

using namespace Dune;
using namespace Fem;

typedef GridSelector::GridType MyGridType;

#ifdef LEAFGRIDPART
// AdaptiveLeafGridPart does not work with YaspGrid for p > 1, since not entities for
// codim > 0 and codim < dim are available, which are needed to build the index set
typedef LeafGridPart< MyGridType > GridPartType;
#else
typedef AdaptiveLeafGridPart< MyGridType > GridPartType;
#endif

typedef FunctionSpace< double, double, GridPartType::dimensionworld, 1 >
FunctionSpaceType;

typedef ExactSolution< FunctionSpaceType > ExactSolutionType;


template <class DiscreteFunctionType>
void writeOut ( VirtualOutStream out, const DiscreteFunctionType &solution )
{
  out << solution;
  out.flush();
}

template <class DiscreteFunctionType>
void readBack ( VirtualInStream in, DiscreteFunctionType &solution )
{
  solution.clear();
  in >> solution;
}

template <class DiscreteFunctionSpaceType>
int algorithm( MyGridType& grid, const int polOrder )
{
  const int step = TestGrid :: refineStepsForHalf();
  grid.globalRefine( 2*step );
  const int rank = grid.comm().rank();

  ////////////////////////////////////////////////////////
  // create data structures (after grid has been refined)
  ////////////////////////////////////////////////////////

  GridPartType gridPart( grid );

  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart, polOrder );

  typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  DiscreteFunctionType solution( "solution", discreteFunctionSpace );
  solution.clear();

  // interpolate
  interpolate( gridFunctionAdapter( ExactSolutionType(), gridPart, discreteFunctionSpace.order() + 2 ), solution );

  // let's check on IO
  DiscreteFunctionType readback( "readback", discreteFunctionSpace );

  std::string casename( "lagrangeinterpolation_r" + std::to_string(rank) + "_p" + std::to_string(polOrder) + "_" );

  BinaryFileOutStream bout( casename + "solution-xdr.tmp" );
  writeOut( virtualize( bout ), solution );

  BinaryFileInStream bin( casename + "solution-xdr.tmp" );
  readBack( virtualize( bin ), readback );
  if( readback != solution )
  {
    std :: cerr << "binary read/write gives different function." << std :: endl;
    assert( readback == solution );
    return 1;
  }

  ASCIIOutStream aout( casename + "solution-ascii.tmp" );
  writeOut( virtualize( aout ), solution );

  ASCIIInStream ain( casename + "solution-ascii.tmp" );
  readBack( virtualize( ain ), readback );
  assert( readback == solution );
  if( readback != solution )
  {
    std :: cerr << "ascii read/write gives different function." << std :: endl;
    assert( readback == solution );
    return 1;
  }

  // output to vtk file
  VTKIO<GridPartType> vtkWriter(gridPart);
  vtkWriter.addVertexData(solution);
  vtkWriter.pwrite(casename+"lagrangeinterpol",
                    Parameter::commonOutputPath().c_str(),"",
                    Dune::VTK::ascii);

  // reset grid
  grid.globalRefine( -2*step );
  return 0;
}

void run( MyGridType& grid, const int polOrder )
{
  // test dynamic Lagrange space
  {
    typedef DynamicLagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType > DiscreteFunctionSpaceType;
    algorithm< DiscreteFunctionSpaceType >(grid, polOrder);
  }

#if HAVE_DUNE_LOCALFUNCTIONS
  // test Lagrange space using LocalFiniteElementSpace
  {
    typedef LagrangeSpace< FunctionSpaceType, GridPartType > DiscreteFunctionSpaceType;
    algorithm< DiscreteFunctionSpaceType >(grid, polOrder);
  }
#endif
}


int main(int argc, char ** argv)
{
  MPIManager :: initialize( argc, argv );
  try
  {
    // add command line parameters to global parameter table
    Dune::Fem::Parameter::append( argc, argv );
    // append parameters from the parameter file
    Dune::Fem::Parameter::append( (argc < 2) ? "parameter" : argv[ 1 ] );

    MyGridType &grid = TestGrid :: grid();

    // if parameter was specified then only run test for this pol order
    // this is primarily for debugging
    if( Dune::Fem::Parameter::exists("fem.lagrange.polynomialOrder") )
    {
      int p = Dune::Fem::Parameter::getValue< int >( "fem.lagrange.polynomialOrder");
      run( grid, p );
    }
    else
    {
      for( int p=1; p <= 3; ++p )
      {
        run( grid, p );
      }
    }

    return 0;
  }
  catch( const Exception& e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}

