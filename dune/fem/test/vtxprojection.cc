#undef NDEBUG

#include <config.h>
#include <iostream>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/operator/projection/vtxprojection.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>

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
    // const int step = TestGrid :: refineStepsForHalf();
    // grid.globalRefine( 2*step );

    GridPartType gridPart( grid );
    DiscreteFunctionSpaceType discreteFunctionSpace( gridPart ); // DG
    typedef LagrangeDiscreteFunctionSpace<FunctionSpaceType, GridPartType,
      1,CachingStorage> LDiscreteFunctionSpaceType;
    LDiscreteFunctionSpaceType lagspace(gridPart);
    ExactSolutionType exactSolution;
    typedef AdaptiveDiscreteFunction<LDiscreteFunctionSpaceType> LagrangeFunctionType;

    // perform the L2Projection
    DiscreteFunctionType solution( "solution", discreteFunctionSpace );
    Fem::L2Projection< ExactSolutionType,  DiscreteFunctionType > dgl2;
    dgl2( exactSolution, solution );

    LagrangeFunctionType contSolution("contSolution",lagspace);
    VtxProjection< DiscreteFunctionType,LagrangeFunctionType > projection;
    projection(solution,contSolution);

    LagrangeFunctionType lagrangeSolution("lagrangeSolution",lagspace);
    Fem::LagrangeInterpolation< ExactSolutionType, LagrangeFunctionType >
      :: interpolateFunction( exactSolution, lagrangeSolution );
    #if 0
    LagrangeType lagrangeContSolution("lagrnageContSolution",lagspace);
    VtxProjection<double,double,LagrangeType,LagrangeType> projection2;
    projection2(lagrangeSolution,lagrangeContSolution);

    // output to vtk file
    VTKIO<GridPartType> vtkWriter(gridPart);
    vtkWriter.addCellData(solution);
    vtkWriter.addVertexData(contSolution);
    lagrangeSolution -= lagrangeContSolution;
    vtkWriter.addVertexData(lagrangeSolution);
    vtkWriter.addVertexData(lagrangeContSolution);
    vtkWriter.pwrite("vtxprojection",
		     Parameter::commonOutputPath().c_str(),".",
		     Dune::VTKOptions::ascii);

#endif
    return 0;
  }
  catch( Exception e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}
