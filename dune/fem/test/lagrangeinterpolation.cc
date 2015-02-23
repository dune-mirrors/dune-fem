#undef NDEBUG
#define ENABLE_ADAPTIVELEAFINDEXSET_FOR_YASPGRID

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
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/io/streams/xdrstreams.hh>
#include <dune/fem/io/streams/asciistreams.hh>
#include <dune/fem/io/streams/virtualstreams.hh>
#include <dune/fem/io/file/vtkio.hh>

#include "testgrid.hh"
#include "exactsolution.hh"
#include "dfspace.hh"

using namespace Dune;
using namespace Fem;

// polynom approximation order of quadratures,
// at least poolynom order of basis functions
#ifdef POLORDER
  const int polOrder = POLORDER;
#else
  const int polOrder = 1;
#endif

typedef GridSelector::GridType MyGridType;

#ifdef LEAFGRIDPART
// AdaptiveLeafGridPart does not work with YaspGrid for p > 1, since not entities for
// codim > 0 and codim < dim are available, which are needed to build the index set
typedef LeafGridPart< MyGridType > GridPartType;
#else
typedef AdaptiveLeafGridPart< MyGridType > GridPartType;
#endif

typedef TestFunctionSpace FunctionSpaceType;

typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
  DiscreteFunctionSpaceType;

typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

typedef ExactSolution< FunctionSpaceType > ExactSolutionType;



void writeOut ( VirtualOutStream out, const DiscreteFunctionType &solution )
{
  out << solution;
  out.flush();
}

void readBack ( VirtualInStream in, DiscreteFunctionType &solution )
{
  solution.clear();
  in >> solution;
}



int main(int argc, char ** argv)
{
  MPIManager :: initialize( argc, argv );
  try
  {
    MyGridType &grid = TestGrid :: grid();
    const int step = TestGrid :: refineStepsForHalf();

    grid.globalRefine( 2*step );

    ////////////////////////////////////////////////////////
    // create data structures (after grid has been refined)
    ////////////////////////////////////////////////////////

    GridPartType gridPart( grid );

    DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
    ExactSolutionType f;
    DiscreteFunctionType solution( "solution", discreteFunctionSpace );
    solution.clear();

    //! perform Lagrange interpolation
    Fem::LagrangeInterpolation< ExactSolutionType, DiscreteFunctionType >
      :: interpolateFunction( f, solution );

    // Let's check on IO
    DiscreteFunctionType readback( "readback", discreteFunctionSpace );

    XDRFileOutStream xout( "solution-xdr.tmp" );
    writeOut( virtualize( xout ), solution );

    XDRFileInStream xin( "solution-xdr.tmp" );
    readBack( virtualize( xin ), readback );
    if( readback != solution )
    {
      std :: cerr << "xdr read/write gives different function." << std :: endl;
      return 1;
    }

    ASCIIOutStream aout( "solution-ascii.tmp" );
    writeOut( virtualize( aout ), solution );

    ASCIIInStream ain( "solution-ascii.tmp" );
    readBack( virtualize( ain ), readback );
    if( readback != solution )
    {
      std :: cerr << "ascii read/write gives different function." << std :: endl;
      return 1;
    }

    // output to vtk file
    VTKIO<GridPartType> vtkWriter(gridPart);
    vtkWriter.addVertexData(solution);
    vtkWriter.pwrite("vtxprojection",
                      Parameter::commonOutputPath().c_str(),"",
                      Dune::VTK::ascii);
    return 0;
  }
  catch( Exception e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}

