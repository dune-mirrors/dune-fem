#undef NDEBUG

#include <config.h>
#include <iostream>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/operator/projection/vtxprojection.hh>
#include <dune/fem/space/common/interpolate.hh>

#if defined USE_BLOCKVECTORFUNCTION
#include <dune/fem/function/blockvectorfunction.hh>
#elif defined USE_VECTORFUNCTION
#include <dune/common/dynvector.hh>
#include <dune/fem/function/vectorfunction.hh>
#else
#include <dune/fem/function/adaptivefunction.hh>
#endif

#include "testgrid.hh"
#include "dfspace.hh"
// #include "vtxl2projection.hh"
#include "exactsolution.hh"
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>

using namespace Dune;
using namespace Fem;

typedef GridSelector::GridType MyGridType;
// typedef HierarchicGridPart< MyGridType > GridPartType;
typedef Fem::AdaptiveLeafGridPart< MyGridType > GridPartType;

typedef TestFunctionSpace FunctionSpaceType;
typedef TestDiscreteFunctionSpace< GridPartType > DiscreteFunctionSpaceType;

#if defined USE_BLOCKVECTORFUNCTION
typedef Fem :: ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >
  DiscreteFunctionType;
#elif defined USE_VECTORFUNCTION
typedef Fem :: ManagedDiscreteFunction
  < Fem::VectorDiscreteFunction
    < DiscreteFunctionSpaceType,
      Dune::DynamicVector< FunctionSpaceType :: RangeFieldType > > >
  DiscreteFunctionType;
#else
typedef Fem :: AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
  DiscreteFunctionType;
#endif

typedef ExactSolution< FunctionSpaceType > ExactSolutionType;

#include <dune/fem/io/file/vtkio.hh>

// main program
int main(int argc, char ** argv)
{
  MPIManager :: initialize( argc, argv );
  try
  {
    MyGridType &grid = TestGrid :: grid();
    const int step = TestGrid :: refineStepsForHalf();
    grid.globalRefine( 2*step );

    for( const auto& entity : elements(grid.leafGridView()) )
    {
      std::remove_reference_t<decltype(entity)>::Geometry::GlobalCoordinate x(0.1);
      x -= entity.geometry().center();
      if( x.two_norm() > 0.3 && x.two_norm() < 0.4 )
        grid.mark( 1, entity );
    }
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();


    GridPartType gridPart( grid );
    DiscreteFunctionSpaceType discreteFunctionSpace( gridPart ); // DG

#if HAVE_DUNE_LOCALFUNCTIONS
    typedef LagrangeSpace<FunctionSpaceType, GridPartType>
      LDiscreteFunctionSpaceType;
    LDiscreteFunctionSpaceType lagspace(gridPart,1);
#else
    typedef LagrangeDiscreteFunctionSpace<FunctionSpaceType, GridPartType, 1>
      LDiscreteFunctionSpaceType;
    LDiscreteFunctionSpaceType lagspace(gridPart);
#endif

    ExactSolutionType exactSolution;
    typedef AdaptiveDiscreteFunction<LDiscreteFunctionSpaceType> LagrangeFunctionType;

    // perform the L2Projection
    DiscreteFunctionType solution( "solution", discreteFunctionSpace );
    interpolate( exactSolution, solution );

    LagrangeFunctionType contSolution("contSolution",lagspace);
    VtxProjection< DiscreteFunctionType,LagrangeFunctionType > projection;
    projection(solution,contSolution);

    LagrangeFunctionType lagrangeSolution("lagrangeSolution",lagspace);
    interpolate( gridFunctionAdapter( exactSolution, gridPart, lagspace.order()+2 ), lagrangeSolution );
#if 0
    LagrangeFunctionType testSolution("test",lagspace);
    VtxProjection<LagrangeFunctionType,LagrangeFunctionType> projection2;
    projection2(lagrangeSolution,testSolution);

    // output to vtk file
    VTKIO<GridPartType> vtkWriter(gridPart);
    vtkWriter.addCellData(solution);
    vtkWriter.addVertexData(contSolution);
    testSolution -= lagrangeSolution;
    vtkWriter.addVertexData(lagrangeSolution);
    vtkWriter.addVertexData(testSolution);
    vtkWriter.pwrite("vtxprojection", Parameter::commonOutputPath().c_str(),".");
#endif
    return 0;
  }
  catch( const Exception& e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}
