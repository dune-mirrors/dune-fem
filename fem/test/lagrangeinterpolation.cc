#include <config.h>
#include <iostream>

#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/common/gridpart.hh>

#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>

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

typedef FunctionSpace< double, double, dimworld, 1 > FunctionSpaceType;

typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
  DiscreteFunctionSpaceType;

typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

typedef ExactSolution< FunctionSpaceType > ExactSolutionType;



int main ()
{
  try
  {
    char tmp[16];
    sprintf( tmp, "%d", dimworld );
    std :: string macroGridName( tmp );
    macroGridName += "dgrid.dgf";

    GridPtr<GridType> gridptr(macroGridName);
    GridType& grid=*gridptr;
    const int step = Dune::DGFGridInfo<GridType>::refineStepsForHalf();
    GridPartType gridPart( grid );

    grid.globalRefine( 2*step );

    DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
    ExactSolutionType f( discreteFunctionSpace );
    DiscreteFunctionType solution( "solution", discreteFunctionSpace );
    solution.clear();

    //! perform Lagrange interpolation
    LagrangeInterpolation< DiscreteFunctionType >
      :: interpolateFunction( f, solution );

    return 0;
  }
  catch( Exception e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}

