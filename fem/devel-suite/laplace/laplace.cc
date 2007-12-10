#include <config.h>

//#define USE_GRAPE HAVE_GRAPE
#ifndef PROBLEM
#warning "PROBLEM not defined. Using SineProblem..."
#define PROBLEM SineProblem
#endif

#ifdef POLORDER
  enum { polynomialOrder = POLORDER };
#else
  enum { polynomialOrder = 1 };
#endif

//- system includes
#include <iostream>

//- dune includes
#include <dune/common/stdstreams.cc>

#include <dune/grid/common/gridpart.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

#if USE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include <dune/fem/space/common/adaptiveleafgridpart.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/misc/l2error.hh>

#include <dune/fem/operator/diffusionoperator.hh>
#include <dune/fem/operator/dirichletboundaryoperator.hh>
#include <dune/fem/operator/cachedlinearoperator.hh>

//- local includes
#include "model.hh"
#include "quadraticproblem.hh"
#include "sineproblem.hh"


using namespace Dune;

typedef double FieldType;

typedef FunctionSpace< FieldType, FieldType, dimworld, 1 > FunctionSpaceType;

typedef PROBLEM< FunctionSpaceType > ProblemType;

typedef ProblemType :: RightHandSideType RightHandSideType;
typedef ProblemType :: ExactSolutionType ExactSolutionType;

typedef LaplaceModel< ProblemType > LaplaceModelType;

typedef LeafGridPart< GridType > GridPartType;

typedef DiscreteFunctionAdapter< RightHandSideType, GridPartType >
  GridRightHandSideType;
typedef DefaultDiffusionOperatorTraits< GridRightHandSideType, polynomialOrder >
  DiffusionOperatorTraitsType;
typedef DiffusionOperator< DiffusionOperatorTraitsType, LaplaceModelType > LaplaceOperatorType;
typedef DirichletBoundaryOperator< LaplaceOperatorType, LaplaceModelType >
  GlobalOperatorType;
typedef CachedLinearOperator< GlobalOperatorType, SparseRowMatrix< FieldType > >
  CachedGlobalOperatorType;

typedef LaplaceOperatorType :: DiscreteFunctionType DiscreteFunctionType;
typedef LaplaceOperatorType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

typedef DiscreteFunctionSpaceType :: GridPartType GridPartType;

typedef OEMCGOp< DiscreteFunctionType, CachedGlobalOperatorType > InverseOperatorType;




FieldType algorithm ( const std :: string &gridFileName, int refinementLevel )
{
  // prepare grid
  GridPtr< GridType > gridPtr( gridFileName );
  gridPtr->globalRefine( refinementLevel );
  
  GridPartType gridPart( *gridPtr );

  // initialize discrete function space
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
  std :: cout << "Solving for " << discreteFunctionSpace.size() << " unknowns."
              << std :: endl << std :: endl;

  // initialize operators
  LaplaceModelType model;

  LaplaceOperatorType laplaceOperator( discreteFunctionSpace, model );
  GlobalOperatorType globalOperator( laplaceOperator, model );
  CachedGlobalOperatorType cachedGlobalOperator( globalOperator );

  double dummy = 0;
  InverseOperatorType inverseOperator( cachedGlobalOperator, dummy, 1e-10, 20000, false );

  // project right hand side
  RightHandSideType rightHandSide( discreteFunctionSpace );
  GridRightHandSideType gridRightHandSide( "continuous right hand side", rightHandSide, gridPart );
  DiscreteFunctionType rhs( "right hand side", discreteFunctionSpace );
  cachedGlobalOperator.rangeProjection()( gridRightHandSide, rhs );

  #if 1
    if( !rhs.dofsValid() )
      std :: cout << "right hand side invalid." << std :: endl;
  #endif

  // solve
  DiscreteFunctionType solution( "solution", discreteFunctionSpace );
  solution.clear();
  inverseOperator( rhs, solution );

  #if USE_GRAPE
    GrapeDataDisplay< GridType > grape( *gridPtr );
    grape.dataDisplay( solution );
  #endif

  ExactSolutionType exactSolution( discreteFunctionSpace );
  L2Error< DiscreteFunctionType > l2error;
  ExactSolutionType :: RangeType error = l2error.norm( exactSolution, solution );
  std :: cout << "L2 error: " << error[ 0 ] << std :: endl << std :: endl;
  
  return error[ 0 ];
}



int main ( int argc, char **argv )
{
  if( argc != 2 )
  {
    std :: cerr << "Usage: " << argv[ 0 ] << " <maxlevel>" << std :: endl;
    return 1;
  }

  int level = atoi( argv[ 1 ] );
  level = (level > 0 ? level - 1 : 0);
  

  std :: string macroGridName( "square.dgf" );

  FieldType error[ 2 ];
  const int steps = DGFGridInfo< GridType > :: refineStepsForHalf();
  for( int i = 0; i < 2; ++i )
    error[ i ] = algorithm( macroGridName, (level+i) * steps );

  const FieldType eoc = log( error[ 0 ] / error[ 1 ] ) / M_LN2;
  std :: cout << "EOC: " << eoc << std :: endl;
  return 0;
}
