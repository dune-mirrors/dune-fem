/**************************************************************************
**       Title: poisson.cc
**    $RCSfile$
**   $Revision$$Name$
**       $Date$
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

// uncomment the following line to use grape
#define USE_GRAPE HAVE_GRAPE

#define VERBOSE false

#include <config.h>

//- system includes
#include <iostream>
#include <sstream>

//- Dune includes 
#include <dune/common/stdstreams.cc>
#include <dune/grid/common/gridpart.hh>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include <dune/fem/space/common/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/solver/inverseoperators.hh>
//#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/misc/l2error.hh>

//- local inlcudes 
#include "laplace.hh"
#include "problem.hh"

#ifndef POLORDER
  #define POLORDER 1
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
typedef RHSFunction< FunctionSpaceType > RHSFunctionType;
typedef ExactSolution< FunctionSpaceType > ExactSolutionType;
typedef Tensor< FunctionSpaceType > TensorType;

//! define the discrete function space our unkown belongs to
typedef LagrangeDiscreteFunctionSpace
  < FunctionSpaceType, GridPartType, POLORDER, CachingStorage >
  DiscreteFunctionSpaceType;

//! define the type of discrete function we are using
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

//! define the discrete laplace operator, see ./fem.cc
typedef LaplaceFEOp< DiscreteFunctionType, TensorType > LaplaceOperatorType;

//! define the inverse operator we are using to solve the system 
// see dune/fem/inverseoperators.hh 
//typedef CGInverseOp < DiscreteFunctionType, LaplaceOperatorType >    InverseOperatorType;
/****************************************/
// or ../../solvers/oemsolver/oemsolvers.hh
typedef OEMCGOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
//typedef OEMBICGSTABOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
//typedef OEMBICGSQOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
//typedef OEMGMRESOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;



//! set the dirichlet points to zero
template< class EntityType, class DiscreteFunctionType >
void boundaryTreatment( const EntityType &entity, DiscreteFunctionType &rhs )
{
  typedef typename DiscreteFunctionType :: FunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;

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
  for( ; it != endit; ++it ) {
    if( !it.boundary() )
      continue;

    LocalFunctionType rhsLocal = rhs.localFunction( entity );
    const LagrangePointSetType &lagrangePointSet
      = discreteFunctionSpace.lagrangePointSet( entity );

    const int face = it.numberInSelf();
    FaceDofIteratorType faceIt
      = lagrangePointSet.template beginSubEntity< faceCodim >( face );
    const FaceDofIteratorType faceEndIt
      = lagrangePointSet.template endSubEntity< faceCodim >( face );
    for( ; faceIt != faceEndIt; ++faceIt )
      rhsLocal[ *faceIt ] = 0;
  }
}



double algorithm ( std :: string &filename, int maxlevel, int turn )
{
  GridPtr< GridType > gridptr( filename ); 
  
  gridptr->globalRefine( maxlevel );

  GridPartType gridPart( *gridptr );

  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
  std::cout << std :: endl << "Solving for " << discreteFunctionSpace.size()
            << " unkowns and polynomial order "
            << DiscreteFunctionSpaceType :: polynomialOrder << "." 
            << std :: endl << std :: endl;
  
  DiscreteFunctionType solution( "solution", discreteFunctionSpace );
  solution.clear();
  DiscreteFunctionType rhs( "rhs", discreteFunctionSpace );
  rhs.clear();
      
  RHSFunctionType f( discreteFunctionSpace ); 
    
  LaplaceOperatorType laplace( discreteFunctionSpace ); 
   
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
    boundaryTreatment( *it , rhs );
  // check the right hand side for NaNs again
  if( !rhs.dofsValid() )
    std :: cout << "right hand side invalid after boundary treatment." << std :: endl;

  // solve the linear system (with CG)
  double dummy = 12345.67890;
  InverseOperatorType cg( laplace, dummy, 1e-8, 20000, VERBOSE );
  cg( rhs, solution );

  // calculation of L2 error
  // polynomial order for this calculation should be higher than the polynomial
  // order of the base functions
  ExactSolutionType u( discreteFunctionSpace ); 
  L2Error< DiscreteFunctionType > l2error;
  DiscreteFunctionSpaceType :: RangeType error = l2error.norm( u, solution );
  std :: cout << "L2 Error: " << error << std :: endl << std :: endl;
   
  #if USE_GRAPE
  // if grape was found then display solution 
  if( turn > 0 )
  {
    GrapeDataDisplay < GridType > grape( *gridptr );
    grape.dataDisplay( solution );
  }
  #endif

  return error[ 0 ];
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
  try {
    if( argc != 2 )
    {
      std :: cerr << "Usage: " << argv[ 0 ] << " <maxlevel>" << std :: endl;
      return 1;
    }
    
    int level = atoi( argv[ 1 ] );
    double error[ 2 ];

    std :: string macroGridName = getMacroGridName( GridType :: dimension );
    std :: cout << "loading dgf: " << macroGridName << std :: endl;
    
    const int step = DGFGridInfo< GridType > :: refineStepsForHalf();
    level = (level > step ? level - step : 0);
    
    for( int i = 0; i < 2; ++i )
      error[ i ] = algorithm( macroGridName, level + i*step, i );

    const double eoc = log( error[ 0 ] / error[ 1 ] ) / M_LN2;
    std :: cout << "EOC = " << eoc << std :: endl;
    
    return 0;
  }
  catch( Exception exception )
  {
    std :: cerr << exception << std :: endl;
    return 1;
  }
}
