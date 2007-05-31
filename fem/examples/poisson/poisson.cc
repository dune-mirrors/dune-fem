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
**************************************************************************/

// uncomment the following line tu use grape
//#define USE_GRAPE HAVE_GRAPE

#define VERBOSE false


//- system includes
#include <iostream>
#include <config.h>

//- Dune includes 
#include <dune/common/stdstreams.cc>
#include <dune/grid/common/gridpart.hh>
#include <dune/grid/io/file/dgfparser/gridtype.hh>

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/operator/inverseoperators.hh>
//#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/misc/l2error.hh>

//- local inlcudes 
#include "laplace.hh"

#ifndef POLORDER
  #define POLORDER 1
#endif



using namespace Dune;

//***********************************************************************
/*! Poisson problem: 

  This is an example how to solve the equation on 
  \f[\Omega = (0,1)^dimworld \f]

  \f[ -\triangle u  = f \ \ \ in \Omega \f]
  \f[  \qquad u = 0  \ \ \ on  \partial\Omega \f]
  \f[ f(x,y,z) = 2 (x-x^2) (y-y^2) +
                 2 (z-z^2) (y-y^2) +    
                 2 (x-x^2) (z-z^2)
  \f]

  An exact solution to this problem is given by 
  \f[ u(x,y,z) = ( x - x^2 ) ( y - y^2 ) ( z - z^2 ) \f]

  with the finite element method using lagrangian elements of polynor order 1.
*/
//***********************************************************************



// forward declaration 
class Tensor; 

//! the index set we are using 
typedef LeafGridPart< GridType > GridPartType;
//typedef LevelGridPart< GridType > GridPartType;

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
typedef FunctionSpace< double, double, dimworld, 1 > FuncSpace;

//! define the discrete function space our unkown belongs to
typedef LagrangeDiscreteFunctionSpace
        < FuncSpace, GridPartType, POLORDER, CachingStorage >
  DiscreteFunctionSpaceType;

//! define the type of discrete function we are using
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

//! define the discrete laplace operator, see ./fem.cc
typedef LaplaceFEOp< DiscreteFunctionType, Tensor > LaplaceOperatorType;

//! define the inverse operator we are using to solve the system 
// see dune/fem/inverseoperators.hh 
//typedef CGInverseOp < DiscreteFunctionType, LaplaceOperatorType >    InverseOperatorType;
/****************************************/
// or ../../solvers/oemsolver/oemsolvers.hh
typedef OEMCGOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
//typedef OEMBICGSTABOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
//typedef OEMBICGSQOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;
//typedef OEMGMRESOp<DiscreteFunctionType,LaplaceOperatorType> InverseOperatorType;



// right hand side of governing problem 
class RHSFunction : public Function< FuncSpace, RHSFunction >
{
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::DomainType DomainType;

private:
  typedef RHSFunction ThisType;
  typedef Function< FuncSpace, ThisType > BaseType;
 
public:
  RHSFunction ( FuncSpace &functionSpace )
  : BaseType( functionSpace )
  {
  }
   
  //  f(x,y,z) = 2 (x-x^2) (y-y^2) +
  //             2 (z-z^2) (y-y^2) +              
  //             2 (x-x^2) (z-z^2)
  void evaluate( const DomainType &x , RangeType &phi ) const
  {
    enum { dimension = DomainType :: dimension };
    
    phi = 0;
    for( int i = 0; i < dimension; ++i ) { 
      RangeType tmp = 2.0;
      for( int j = 1; j < dimension; ++j ) {
        const int idx = (i + j) % dimension;
        tmp *= x[ idx ] - SQR( x[ idx ] );
      }
      phi += tmp;
    }
  }
};



//! the exact solution to the problem for EOC calculation 
class ExactSolution : public Function < FuncSpace , ExactSolution > 
{
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::RangeFieldType RangeFieldType;
  typedef FuncSpace::DomainType DomainType;
public:
  ExactSolution (FuncSpace &f) : Function < FuncSpace , ExactSolution > ( f ) {}
 
  //! u(x,y,z) = (x-x^2)*(y-y^2)*(z-z^2)
  void evaluate (const DomainType & x , RangeType & ret) const
  {
    ret = 1.0;
    for(int i=0; i<DomainType::dimension; i++)
      ret *= ( x[i] - SQR(x[i]) );
  }
  void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
  {
    evaluate ( x , ret );
  }
};



// diffusion coefficient for this problem the id 
class Tensor : public Function < FunctionSpace < double , double, dimworld , 1 >, Tensor >
{
  typedef FunctionSpace< double , double, dimworld , 1 >  FuncSpace;
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::DomainType DomainType;

public:
  // Constructor 
  Tensor (FuncSpace &f)
    : Function < FuncSpace , Tensor > ( f )  { } ;
  // eval Tensor 
  void evaluate (int i, int j, const DomainType & x1 , RangeType & ret) const
  {
    evaluate(x1,0.0,ret);
  }
  void evaluate (const DomainType & x1 , RangeType & ret) const
  {
    evaluate(x1,0.0,ret);
  }
  void evaluate (const DomainType & x1 ,double time, RangeType & ret) const
  {
    ret[0] = 1.0;
  }
};//end class Tensor



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
      
  RHSFunction f( discreteFunctionSpace ); 
    
  LaplaceOperatorType laplace( discreteFunctionSpace, 
                               LaplaceOperatorType :: ASSEMBLED );
   
   //! build right hand side, does not allocate b!
  RightHandSideAssembler< DiscreteFunctionType >
    :: assemble< 2 * DiscreteFunctionSpaceType :: polynomialOrder >( f , rhs );
    
  // set Dirichlet Boundary to zero 
  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType; 
  IteratorType endit = discreteFunctionSpace.end();
  for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
    boundaryTreatment( *it , rhs );

  //laplace.print();
  //rhs.print( std::cout );
   
  // solve the linear system (with CG)
  double dummy = 12345.67890;
  InverseOperatorType cg( laplace, dummy, 1e-8, 20000, VERBOSE );
  cg( rhs, solution );

  // calculation of L2 error
  // polynomial order for this calculation should be higher than the polynomial
  // order of the base functions
  ExactSolution u( discreteFunctionSpace ); 
  L2Error< DiscreteFunctionType > l2error;
  DiscreteFunctionSpaceType :: RangeType error = l2error.norm( u, solution );
  std :: cout << "L2 Error: " << error << std :: endl << std :: endl;
   
  #if USE_GRAPE
  // if grape was found then display solution 
  if( turn > 0 ) {
    GrapeDataDisplay < GridType > grape( *gridptr );
    grape.dataDisplay( solution );
  }
  #endif

  return error[ 0 ];
}



// main programm, run algorithm twice to calc EOC 
int main( int argc, char **argv )
{
  if( argc != 2 ) {
    std :: cerr << "Usage: " << argv[ 0 ] << " <maxlevel>" << std :: endl;
    return 1;
  }
  
  int level = atoi( argv[ 1 ] );
  double error[ 2 ];

  std :: string macroGridName( "square.dgf" );
  std :: cout << "loading dgf: " << macroGridName << std :: endl;
  
  const int step = DGFGridInfo< GridType > :: refineStepsForHalf();
  level = (level > step ? level - step : 0);
  
  for( int i = 0; i < 2; ++i )
    error[ i ] = algorithm( macroGridName, level + i*step, i );

  const double eoc = log( error[ 0 ] / error[ 1 ] ) / M_LN2;
  std :: cout << "EOC = " << eoc << std :: endl;
  
  return 0;
}
