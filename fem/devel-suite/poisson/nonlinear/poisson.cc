/**************************************************************************
**       Title: poisson.cc
**    $RCSfile$
**   $Revision: 3976 $$Name$
**       $Date: 2008-09-16 15:24:56 +0200 (Tue, 16 Sep 2008) $
**   Copyright: GPL Author: robertk 
** Description: File demonstrating a simple numerics problem on arbitrary
**              grids: poisson-problem with known solution is given
**              and compared with numerical solution (EOC)
**              Dune grid parser is used.
**
**              For changing grid-types, compile with
**
**              make clean
**
**              and one of the following
**           
**              make GRIDTYPE=YASPGRID       (default)
**                    -> compiles and works correctly
**              make
**                    -> compiles and works correctly
**              make GRIDTYPE=ALBERTAGRID
**                    -> compiles and works correctly
**              make GRIDTYPE=SGRID
**                    -> compiles and works correctly
**              make GRIDTYPE=ALUGRID_SIMPLEX
**                    -> compiles and works correctly
**
**
**              Similarly, the polynomial order can be specified by
**
**              make POLORDER=2
**
**************************************************************************/

#define VERBOSE false

#define USE_TWISTFREE_MAPPER

#include <config.h>

#ifdef ENABLE_TIMING
#include <time.h>
#endif

//- system includes
#include <iostream>
#include <sstream>

//- Dune includes 
#include <dune/common/stdstreams.cc>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/operator/matrix/blockmatrix.hh>
#include <dune/fem/operator/matrix/ontheflymatrix.hh>
#include <dune/fem/operator/matrix/istlmatrix.hh>
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/fem/solver/istlsolver.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/io/visual/grape/datadisp/errordisplay.hh>

#include <dune/fem/misc/mpimanager.hh>

//- local inlcudes 
#include "laplace.hh"
#if (PROBLEM==2)
#include "cornerproblem.hh"
#elif (PROBLEM==1)
#include "albertaproblem.hh"
#else
#include "problem.hh"
#endif

#ifndef POLORDER
  #define POLORDER 1
#endif

#ifndef USE_GRAPE
#define USE_GRAPE 0
#endif

//***********************************************************************
/*! Poisson problem: 

  This is an example solving the poisson problem
  \f{displaymath}
  \begin{array}{rcll}
  -\triangle u &=& f & \quad \mbox{in}\ \Omega\\
             u &=& 0 & \quad \mbox{on}\ \partial\Omega
  \end{array}
  \f}
  with the finite element method using Lagrangian elements. The polynomial
  order is given by POLORDER.
  
  In this example, $\Omega = ]0,1[^{dimworld}$ and
  \f[
  f( x, y, z ) = 2 \sum_{i=1}^{dimworld} \prod_{j \neq i} (x_j-x_j^2)
  \f]

  The exact solution to the poisson problem is then given by
  \f[
  u( x, y, z ) = \prod_{i=1}^{dimworld} (x_i - x_i^2).
  \f]
*/
//***********************************************************************

using namespace Dune;



/** choose the grid partition (and hence the index set) to use
 *
 *  \note Not all index sets a continuous. The LeafIndexSet for AlbertaGrid,
 *        for example, is not. If you want to use OEM solvers, the index set
 *        must be continuous. In such a case use AdaptiveLeafGridPart.
 */
//typedef LeafGridPart< GridType > GridPartType;
//typedef LevelGridPart< GridType > GridPartType;
typedef AdaptiveLeafGridPart< GridType > GridPartType;

//! define the function space, \f[ \R^n \rightarrow \R \f]
// see dune/common/functionspace.hh
typedef FunctionSpace< double, double, dimworld, 1 > FunctionSpaceType;

// The data functions (as defined in problem.hh)
typedef MassFunction< FunctionSpaceType >MassFunctionType;
typedef ExactSolution< FunctionSpaceType > ExactSolutionType;
typedef RHSFunction< FunctionSpaceType, ExactSolutionType, MassFunctionType> RHSFunctionType;
//typedef ExactSolution< FunctionSpaceType > ExactSolutionType;
typedef Tensor< FunctionSpaceType > TensorType;

typedef DiscreteFunctionAdapter< ExactSolutionType, GridPartType >
  GridExactSolutionType;

//! define the discrete function space our unkown belongs to
typedef LagrangeDiscreteFunctionSpace
  < FunctionSpaceType, GridPartType, POLORDER, CachingStorage >
  DiscreteFunctionSpaceType;

//! define the type of discrete function we are using
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
//typedef BlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
//typedef ManagedDiscreteFunction< VectorDiscreteFunction< DiscreteFunctionSpaceType, DynamicVector< double > > > DiscreteFunctionType;

//! define the type of the system matrix object
typedef SparseRowMatrixTraits < DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > MatrixObjectTraits;
//typedef BlockMatrixTraits < DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > MatrixObjectTraits;
//typedef OnTheFlyMatrixTraits< DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > MatrixObjectTraits;
//typedef ISTLMatrixTraits < DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > MatrixObjectTraits;

//! define the discrete laplace operator, see ./fem.cc
typedef LaplaceFEOp< DiscreteFunctionType, MatrixObjectTraits, TensorType,MassFunctionType >
  LaplaceOperatorType;

// typedef LaplaceFEOp< DiscreteFunctionType, MatrixObjectTraits, TensorType >
//   LaplaceOperatorType;

//! define the inverse operator we are using to solve the system 
typedef CGInverseOp< DiscreteFunctionType, LaplaceOperatorType >
//typedef ISTLBICGSTABOp< DiscreteFunctionType, LaplaceOperatorType >
//typedef OEMCGOp<DiscreteFunctionType,LaplaceOperatorType>
//typedef OEMBICGSTABOp<DiscreteFunctionType,LaplaceOperatorType>
//typedef OEMBICGSQOp<DiscreteFunctionType,LaplaceOperatorType>
//typedef OEMGMRESOp<DiscreteFunctionType,LaplaceOperatorType>
  InverseOperatorType;



//! set the dirichlet points to zero
template< class EntityType, class GridFunctionType, class DiscreteFunctionType >
void boundaryTreatment( const EntityType &entity,
                        const GridFunctionType &exactSolution,
                        DiscreteFunctionType &rhs )
{
  typedef typename DiscreteFunctionType :: FunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

  typedef typename GridFunctionType :: LocalFunctionType LocalExactSolutionType;

  typedef typename DiscreteFunctionSpaceType :: LagrangePointSetType
    LagrangePointSetType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  
  enum { faceCodim = 1 };
  typedef typename GridPartType :: IntersectionIteratorType
    IntersectionIteratorType;
  typedef typename LagrangePointSetType :: template Codim< faceCodim > 
                                        :: SubEntityIteratorType
    FaceDofIteratorType;

  const DiscreteFunctionSpaceType &discreteFunctionSpace = rhs.space();
  const GridPartType &gridPart = discreteFunctionSpace.gridPart();
    
  IntersectionIteratorType it = gridPart.ibegin( entity );
  const IntersectionIteratorType endit = gridPart.iend( entity );
  for( ; it != endit; ++it )
  {
    if( !it->boundary() )
      continue;

    LocalExactSolutionType exactLocal = exactSolution.localFunction( entity );
    LocalFunctionType rhsLocal = rhs.localFunction( entity );
    const LagrangePointSetType &lagrangePointSet
      = discreteFunctionSpace.lagrangePointSet( entity );

    const int face = it->numberInSelf();
    FaceDofIteratorType faceIt
      = lagrangePointSet.template beginSubEntity< faceCodim >( face );
    const FaceDofIteratorType faceEndIt
      = lagrangePointSet.template endSubEntity< faceCodim >( face );
    for( ; faceIt != faceEndIt; ++faceIt )
    {
      typename LocalExactSolutionType :: RangeType phi;
      exactLocal.evaluate( lagrangePointSet.point( *faceIt ), phi );
      rhsLocal[ *faceIt ] = phi[ 0 ];
    }
  }
}



void solve ( LaplaceOperatorType &laplace,
             const DiscreteFunctionType &rhs,
             DiscreteFunctionType &solution )
{
#ifdef ENABLE_TIMING
  time_t starttime = time( NULL );
#endif
  
  // solve the linear system (with CG)
  double dummy = 12345.67890;
  InverseOperatorType cg( laplace, dummy, 1e-8, 20000, VERBOSE );
  cg( rhs, solution );

#ifdef ENABLE_TIMING
  time_t endtime = time( NULL );
  std :: cout << "Time needed by solver: " << (endtime - starttime) << std :: endl;
#endif
}



double algorithm ( GridType &grid, int turn )
{
  GridPartType gridPart( grid );

  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
  std::cout << std :: endl << "Solving for " << discreteFunctionSpace.size()
            << " unkowns and polynomial order "
            << DiscreteFunctionSpaceType :: polynomialOrder << "." 
            << std :: endl << std :: endl;

  ExactSolutionType u( discreteFunctionSpace ); 

  GridExactSolutionType ugrid( "exact solution", u, gridPart, DiscreteFunctionSpaceType :: polynomialOrder + 1 );

  DiscreteFunctionType solution( "solution", discreteFunctionSpace );
  solution.clear();
  DiscreteFunctionType rhs( "rhs", discreteFunctionSpace );
  rhs.clear();

  MassFunctionType mass( discreteFunctionSpace);
  RHSFunctionType f( discreteFunctionSpace, u, mass ); 
    
  LaplaceOperatorType laplace( discreteFunctionSpace, mass );   
//   LaplaceOperatorType laplace( discreteFunctionSpace ); 
   
   //! build right hand side, does not allocate b!
  RightHandSideAssembler< DiscreteFunctionType >
    :: assemble< 2 * DiscreteFunctionSpaceType :: polynomialOrder >( f , rhs );
  // we're having some trouble with NaNs lately, so check the right hand side for NaNs
  if( !rhs.dofsValid() )
    std :: cout << "right hand side invalid before boundary treatment." << std :: endl;
   
  // set Dirichlet Boundary to zero 
  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType; 
  IteratorType endit = discreteFunctionSpace.end();
  for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
  {
    boundaryTreatment( *it, ugrid, rhs );
    boundaryTreatment( *it, ugrid, solution );
  }
  // check the right hand side for NaNs again
  if( !rhs.dofsValid() )
    std :: cout << "right hand side invalid after boundary treatment." << std :: endl;

  solve( laplace, rhs, solution );

#if 0
  L2Norm< GridPartType > norm( gridPart );
  double error = norm.distance( ugrid, solution );
  std :: cout << "L2-Error: " << error << std :: endl << std :: endl;
#else
  H1Norm< GridPartType > norm( gridPart );
  double error = norm.distance( ugrid, solution );
  std :: cout << "H1-Error: " << error << std :: endl << std :: endl;
#endif

#if (USE_GRAPE && HAVE_GRAPE)
  // if grape was found then display solution
  if( turn > 0 )
  {
    GrapeDataDisplay < GridType > grape( grid );
    DisplayErrorFunction< DiscreteFunctionType, ExactSolutionType, false>
      displayErrorFunction( grape, solution, u );
    grape.addData( solution );
    grape.display();
  }
#endif

  //return error[ 0 ];
  return error;
}



double algorithm ( std :: string &filename, int maxlevel, int turn )
{
  GridPtr< GridType > gridptr( filename ); 
  
  gridptr->loadBalance();

  gridptr->globalRefine( maxlevel );

  return algorithm( *gridptr, turn );
}



std :: string getMacroGridName( unsigned int dimension )
{
  std :: ostringstream s;
  s << dimension << "dgrid.dgf";
  return s.str();
}



// main programm, run algorithm twice to calc EOC 
int main( int argc, char **argv )
{
  // initialize MPI 
  MPIManager :: initialize( argc, argv );
  const int rank = MPIManager :: rank ();
  
  try
  {
    if( argc < 2 )
    {
      if( rank == 0 )
        std :: cerr << "Usage: " << argv[ 0 ] << " <maxlevel> [macrogrid]"
                    << std :: endl;
      return 1;
    }
    
    int level = atoi( argv[ 1 ] );
    double error[ 2 ];

    std :: string macroGridName
      = (argc > 2 ? argv[ 2 ] : getMacroGridName( GridType :: dimension ));
    if( rank == 0 )
      std :: cout << "loading macro grid: " << macroGridName << std :: endl;
    
    GridPtr< GridType > gridptr( macroGridName ); 
  
    gridptr->loadBalance();

    const int step = DGFGridInfo< GridType > :: refineStepsForHalf();
    level = (level > step ? level - step : 0);
    
    gridptr->globalRefine( level * step );
    
    for( int i = 0; i < 2; ++i )
    {
      gridptr->globalRefine( step * i );
      error[ i ] = algorithm( *gridptr, i );
    }

    const double eoc = log( error[ 0 ] / error[ 1 ] ) / M_LN2;
    std :: cout << "EOC = " << eoc << std :: endl;
    
    return 0;
  }
  catch( Exception &exception )
  {
    if( rank == 0 )
      std :: cerr << exception << std :: endl;
    return 1;
  }
}
