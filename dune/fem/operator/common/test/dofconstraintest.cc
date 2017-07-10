#include <config.h>

#include <dune/grid/yaspgrid.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/io/file/vtkio.hh>

#include <dune/fem/test/testgrid.hh>
#include <dune/fem/test/exactsolution.hh>

#include <dune/fem/operator/common/dofcontraints.hh>

typedef Dune::GridSelector::GridType HGridType;
typedef Dune::Fem::AdaptiveLeafGridPart< HGridType > GridPartType;
typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimension, 2 > FunctionSpaceType;
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 3 > DiscreteFunctionSpaceType;
typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
typedef Dune::Fem::ExactSolution< FunctionSpaceType > ExactSolutionType;

int main(int argc, char ** argv)
{
  Dune::Fem::MPIManager :: initialize( argc, argv );
  // try
  {
    auto &grid = Dune::Fem::TestGrid :: grid();
    const int step = Dune::Fem::TestGrid :: refineStepsForHalf();
    grid.globalRefine( 2*step );
    GridPartType gridPart( grid );
    DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
    DiscreteFunctionType solution( "solution", discreteFunctionSpace );

    solution.clear();
    interpolate( gridFunctionAdapter( ExactSolutionType(), gridPart, discreteFunctionSpace.order() + 2 ), solution );

    ConstrainOnBoundary< DiscreteFunctionSpaceType > mask( discreteFunctionSpace );
    DofConstraints< DiscreteFunctionSpaceType, ConstrainOnBoundary< DiscreteFunctionSpaceType >  >
        constrain( discreteFunctionSpace, mask );
    constrain( solution );

    Dune::Fem::VTKIO<GridPartType> vtkWriter(gridPart);
    vtkWriter.addVertexData(solution);
    vtkWriter.pwrite("constrain", Dune::Fem::Parameter::commonOutputPath().c_str(),".");
    return 0;
  }
#if 0
  catch( Exception e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
#endif
}
