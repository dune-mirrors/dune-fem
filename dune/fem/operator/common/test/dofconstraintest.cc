#include <config.h>

#include <dune/grid/yaspgrid.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/localfunctionadapter.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/io/file/vtkio.hh>

#include <dune/fem/test/testgrid.hh>
#include <dune/fem/test/exactsolution.hh>

#include <dune/fem/operator/common/dofcontraints.hh>

#include <dune/fem/operator/common/test/model.hh>

typedef Dune::GridSelector::GridType HGridType;
typedef Dune::Fem::AdaptiveLeafGridPart< HGridType > GridPartType;
typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimension, 1 > FunctionSpaceType;
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
    interpolate( gridFunctionAdapter( ExactSolutionType(), gridPart, discreteFunctionSpace.order() + 2 ), solution );

    DiscreteFunctionType interpol( "interpolation", discreteFunctionSpace );
    interpol.assign( solution );

    Dune::Fem::ConstrainOnFullBoundary< DiscreteFunctionSpaceType > mask( discreteFunctionSpace );
    Dune::Fem::DofConstraints< DiscreteFunctionSpaceType, Dune::Fem::ConstrainOnFullBoundary< DiscreteFunctionSpaceType >  >
        constrain( discreteFunctionSpace, mask );
    constrain( solution );
    DiscreteFunctionType zerobc( "zerobc", discreteFunctionSpace );
    zerobc.assign( solution );

    constrain.set( gridFunctionAdapter( ExactSolutionType(), gridPart, discreteFunctionSpace.order() + 2 ) );
    constrain( solution );
    solution -= interpol;

    DiscreteFunctionType modelTest( "modelTest", discreteFunctionSpace );
    interpolate( gridFunctionAdapter( ExactSolutionType(), gridPart, discreteFunctionSpace.order() + 2 ), modelTest );
    typedef Model< FunctionSpaceType, GridPartType > ModelType;
    ModelType model;
    Dune::Fem::ConstrainOnBoundary< DiscreteFunctionSpaceType, ModelType > modelMask( discreteFunctionSpace, model );
    Dune::Fem::DofConstraints< DiscreteFunctionSpaceType, Dune::Fem::ConstrainOnBoundary< DiscreteFunctionSpaceType, ModelType >  >
        modelConstrain( discreteFunctionSpace, modelMask );
    Dune::Fem::PiecewiseGridFunction<ModelType> lpw( model, gridPart, 5 );
    Dune::Fem::LocalFunctionAdapter < Dune::Fem::PiecewiseGridFunction<ModelType> >
      pwgf( "pwgf", lpw, gridPart, 5 );
    modelConstrain.set( pwgf );
    modelConstrain( modelTest );

    Dune::Fem::VTKIO<GridPartType> vtkWriter(gridPart);
    vtkWriter.addVertexData(solution);
    vtkWriter.addVertexData(zerobc);
    vtkWriter.addVertexData(interpol);
    vtkWriter.addVertexData(modelTest);
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
