#undef NDEBUG

#include <config.h>
#include <iostream>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/operator/projection/vtxprojection.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>

#if defined USE_BLOCKVECTORFUNCTION
#include <dune/fem/function/blockvectorfunction.hh>
#elif defined USE_VECTORFUNCTION
#include <dune/fem/storage/vector.hh>
#include <dune/fem/function/vectorfunction.hh>
#elif defined USE_ATTACHEDFUNCTION
#include <dune/fem/function/attachedfunction/function.hh>
#else
#include <dune/fem/function/adaptivefunction.hh>
#endif

#include "testgrid.hh"
#include "dfspace.hh"
#include "dgl2projection.hh"
// #include "vtxl2projection.hh"
#include "exactsolution.hh"

using namespace Dune;
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
// typedef HierarchicGridPart< GridType > GridPartType;
typedef AdaptiveLeafGridPart< GridType > GridPartType;

typedef TestFunctionSpace FunctionSpaceType;
typedef TestDiscreteFunctionSpace< GridPartType > DiscreteFunctionSpaceType;

#if defined USE_BLOCKVECTORFUNCTION
typedef BlockVectorDiscreteFunction< DiscreteFunctionSpaceType >
  DiscreteFunctionType;
#elif defined USE_VECTORFUNCTION
typedef ManagedDiscreteFunction
  < VectorDiscreteFunction
    < DiscreteFunctionSpaceType,
      DynamicVector< FunctionSpaceType :: RangeFieldType > > >
  DiscreteFunctionType;
#elif defined USE_ATTACHEDFUNCTION
typedef AttachedDiscreteFunction< DiscreteFunctionSpaceType >
  DiscreteFunctionType;
#else
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
  DiscreteFunctionType;
#endif

typedef ExactSolution< FunctionSpaceType > ExactSolutionType;

#include <dune/fem/io/file/vtkio.hh>

int main(int argc, char ** argv) 
{
  MPIManager :: initialize( argc, argv );
  try
  {
    GridType &grid = TestGrid :: grid();
    // const int step = TestGrid :: refineStepsForHalf();
    // grid.globalRefine( 2*step );

    GridPartType gridPart( grid );
    DiscreteFunctionSpaceType discreteFunctionSpace( gridPart ); // DG
    typedef LagrangeDiscreteFunctionSpace<FunctionSpaceType, GridPartType,
      1,CachingStorage> LDiscreteFunctionSpaceType; 
    LDiscreteFunctionSpaceType lagspace(gridPart);
    ExactSolutionType exactSolution( lagspace );
    typedef AdaptiveDiscreteFunction<LDiscreteFunctionSpaceType> LagrangeFunctionType;

    // perform the L2Projection
    DiscreteFunctionType solution( "solution", discreteFunctionSpace );
    DGL2Projection< DiscreteFunctionType > :: project( exactSolution, solution );
    
    LagrangeFunctionType contSolution("contSolution",lagspace);
    VtxProjection<double,double,DiscreteFunctionType,LagrangeFunctionType> projection;
    projection(solution,contSolution);

    LagrangeFunctionType lagrangeSolution("lagrangeSolution",lagspace);
    LagrangeInterpolation< LagrangeFunctionType >
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
