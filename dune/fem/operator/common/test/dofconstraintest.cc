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

#include <dune/fem/operator/common/dofconstraints.hh>

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

    // first simply interpolate the given function (store in interpol)
    DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
    DiscreteFunctionType solution( "solution", discreteFunctionSpace );
    interpolate( gridFunctionAdapter( ExactSolutionType(), gridPart, discreteFunctionSpace.order() + 2 ), solution );
    DiscreteFunctionType interpol( "interpolation", discreteFunctionSpace );
    interpol.assign( solution );

    // next use zero b.c. on the whole boundary and store in zerobc
    Dune::Fem::ConstrainOnFullBoundary< DiscreteFunctionSpaceType > mask( discreteFunctionSpace );
    Dune::Fem::DofConstraints< DiscreteFunctionSpaceType, Dune::Fem::ConstrainOnFullBoundary< DiscreteFunctionSpaceType >  >
        constraints( discreteFunctionSpace, mask );
    constraints.constrain( solution );
    DiscreteFunctionType zerobc( "zerobc", discreteFunctionSpace );
    zerobc.assign( solution );

    // now use the value from the given function as b.c. on the whole
    // domain then solution == interpol (so after this solution==0)
    constraints.set( gridFunctionAdapter( ExactSolutionType(), gridPart, discreteFunctionSpace.order() + 2 ) );
    constraints.constrain( solution );
    solution -= interpol;

    // now test a b.c. given by the 'model'
    DiscreteFunctionType modelTest( "modelTest", discreteFunctionSpace );
    interpolate( gridFunctionAdapter( ExactSolutionType(), gridPart, discreteFunctionSpace.order() + 2 ), modelTest );

    typedef Model< FunctionSpaceType, GridPartType > ModelType;
    typedef Dune::Fem::ConstrainOnBoundary< DiscreteFunctionSpaceType, ModelType > ModelBoundaryMask;
    typedef Dune::Fem::DofConstraints< DiscreteFunctionSpaceType, ModelBoundaryMask > DirichletConstraints;
    ModelType model;
    ModelBoundaryMask  modelMask( discreteFunctionSpace, model );
    DirichletConstraints modelConstraints( discreteFunctionSpace, modelMask );
    Dune::Fem::PiecewiseGridFunction<ModelType> lpw( model, gridPart, 5 );
    Dune::Fem::LocalFunctionAdapter < Dune::Fem::PiecewiseGridFunction<ModelType> >
      pwgf( "pwgf", lpw, gridPart, 5 );
    modelConstraints.set( pwgf );
    modelConstraints.constrain( modelTest );

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
