#include <config.h>
#include <iostream>

#include <dune/grid/common/gridpart.hh>

#include <dune/fem/misc/double.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/io/streams/xdrstreams.hh>
#include <dune/fem/io/streams/asciistreams.hh>

#include "testgrid.hh"
#include "exactsolution.hh"

using namespace Dune;

// polynom approximation order of quadratures, 
// at least poolynom order of basis functions 
#ifdef POLORDER
  const int polOrder = POLORDER;
#else
  const int polOrder = 1;
#endif

typedef LeafGridPart< GridType > GridPartType;

typedef FunctionSpace< double, Double, dimworld, 1 > FunctionSpaceType;

typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
  DiscreteFunctionSpaceType;

typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

typedef ExactSolution< FunctionSpaceType > ExactSolutionType;



int main ()
{
  try
  {
    GridType &grid = TestGrid :: grid();
    const int step = TestGrid :: refineStepsForHalf();
    GridPartType gridPart( grid );

    grid.globalRefine( 2*step );

    DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
    ExactSolutionType f( discreteFunctionSpace );
    DiscreteFunctionType solution( "solution", discreteFunctionSpace );
    solution.clear();

    std :: cout << "maxDofs = " << discreteFunctionSpace.mapper().numDofs() << std :: endl;

    //! perform Lagrange interpolation
    LagrangeInterpolation< DiscreteFunctionType >
      :: interpolateFunction( f, solution );

    // Let's check on IO
    DiscreteFunctionType readback( "readback", discreteFunctionSpace );
    
    XDRFileOutStream out( "solution-xdr.tmp" );
    out << solution;
    out.flush();
    readback.clear();
    XDRFileInStream in( "solution-xdr.tmp" );
    in >> readback;
    if( readback != solution )
      return 1;

    ASCIIOutStream aout( "solution-ascii.tmp" );
    aout << solution;
    aout.flush();
    readback.clear();
    ASCIIInStream ain( "solution-ascii.tmp" );
    ain >> readback;
    if( readback != solution )
      return 1;

    return 0;
  }
  catch( Exception e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}

