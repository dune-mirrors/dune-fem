#include <config.h>

// to use grape, set to WANT_GRAPE to 1
#define WANT_GRAPE 1
#define SHOW_INTERPOLATION 0
#define SHOW_RESTRICT_PROLONG 1

// to write out the data, set WRITE_DATA to 1
#define WRITE_DATA 0

// to use generic adaption, set GENERIC_ADAPT to 1
#define GENERIC_ADAPT 1

// polynomial order of base functions
const int polOrder = POLORDER;

#include <iostream>
#include <dune/common/stdstreams.cc>

#include <dune/grid/common/gridpart.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

#include <dune/fem/space/common/adaptiveleafgridpart.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/space/lagrangespace/adaptmanager.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachequad.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>

#ifndef GRIDDIM 
  #define GRIDDIM dimworld 
#endif

#if HAVE_GRAPE
//#if HAVE_GRAPE && (GRIDDIM > 1)
  #define USE_GRAPE WANT_GRAPE
#else
  #define USE_GRAPE 0
  #if WANT_GRAPE
    #warning "Grape was not found by configure."
  #endif
#endif

#if USE_GRAPE 
  #include <dune/grid/io/visual/grapedatadisplay.hh>
#endif
#include <dune/fem/io/file/grapedataio.hh>


using namespace Dune;

//! type of the grid partition we are using 
typedef AdaptiveLeafGridPart< GridType > GridPartType;
// typedef HierarchicGridPart< GridType > GridPartType;

//! type of the function space
typedef FunctionSpace< double, double, GRIDDIM, 1 > FunctionSpaceType;

//! type of the discrete function space our unkown belongs to
typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
  DiscreteFunctionSpaceType;

//! type of the discrete function we are using
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;



//! type of the DoF manager
typedef DofManager< GridType > DofManagerType;
//! type of the DoF manager factory
typedef DofManagerFactory< DofManagerType > DofManagerFactoryType;

//! type of restrict-prolong operator
typedef RestrictProlongDefault< DiscreteFunctionType >
  RestrictProlongOperatorType;
//! type of the adaption manager
typedef AdaptationManager< GridType, RestrictProlongOperatorType >
  AdaptationManagerType;



//! exact function we will interpolate
class ExactSolution
: public Function< FunctionSpaceType, ExactSolution >
{
public:
  typedef FunctionSpaceType :: DomainFieldType DomainFieldType;
  typedef FunctionSpaceType :: RangeFieldType RangeFieldType;

  typedef FunctionSpaceType :: DomainType DomainType;
  typedef FunctionSpaceType :: RangeType RangeType;
  typedef FunctionSpaceType :: JacobianRangeType JacobianRangeType;

private:
  typedef ExactSolution ThisType;
  typedef Function< FunctionSpaceType, ThisType > BaseType;

public:
  ExactSolution( FunctionSpaceType &functionSpace )
  : BaseType( functionSpace )
  {
  }

  void evaluate( const DomainType &x, RangeType &phi ) const
  {
    phi = 1;
    for( int i = 0; i < DomainType :: dimension; ++i )
      // phi[ 0 ] += x[ i ] * x[ i ]; 
      phi[ 0 ] *= sin( M_PI * x[ i ] ); 
  }

  void evaluate( const DomainType &x, RangeFieldType t, RangeType &phi ) const
  {
    evaluate( x, phi );
  }

  void jacobian( const DomainType &x, JacobianRangeType &Dphi ) const
  {
    Dphi = 1;
    for( int i = 0; i < DomainType :: dimension; ++i )
      for( int j = 0; j < DomainType :: dimension; ++j )
        // Dphi[ 0 ][ j ] *= ((i != j) ? 1. : 2.*x[i]);
        Dphi[ 0 ][ j ] *= ((i != j) ? sin( M_PI * x[ i ]) : M_PI * cos( M_PI * x[ i ] ));
  }

  void jacobian( const DomainType &x, RangeFieldType t, JacobianRangeType &Dphi ) const
  {
    jacobian( x, Dphi );
  }
};
 


// calculates \Vert u - u_h \Vert_{L^2}
template< class DiscreteFunctionImp >
class L2Error
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
    
  typedef typename DiscreteFunctionType :: FunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType :: LocalFunctionType
    LocalFunctionType;
  
  typedef typename DiscreteFunctionSpaceType :: DomainFieldType
    DomainFieldType;
  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType
    RangeFieldType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  
  enum { DimRange = DiscreteFunctionSpaceType :: DimRange };

  typedef typename GridPartType :: GridType GridType;
  typedef typename GridType :: template Codim< 0 > :: Entity Entity0Type;
  typedef typename Entity0Type :: Geometry GeometryType;
 
  
public:
  template< class FunctionType >
  static RangeFieldType norm ( const FunctionType &function,
                               DiscreteFunctionType &discreteFunction,
                               double time,
                               int polOrder )
  {
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

    const DiscreteFunctionSpaceType &discreteFunctionSpace
      = discreteFunction.space();

    RangeFieldType error( 0 );
    IteratorType endit = discreteFunctionSpace.end();
    for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
    {
      CachingQuadrature< GridPartType, 0 > quadrature( *it, polOrder );
      LocalFunctionType localFunction = discreteFunction.localFunction( *it );

      const GeometryType &geometry = (*it).geometry();
      
      const int numQuadraturePoints = quadrature.nop();
      for( int qp = 0; qp < numQuadraturePoints; ++qp )
      {
        const DomainType &x = quadrature.point( qp );
        const double weight
          = quadrature.weight( qp ) * geometry.integrationElement( x );

        RangeType phi, psi;
        function.evaluate( geometry.global( x ), time, phi );
        localFunction.evaluate( quadrature[ qp ], psi );

        for( int i = 0; i < DimRange; ++i )
          error += weight * SQR( phi[ i ] - psi[ i ] );
      }
    }
    
    return sqrt( error );
  }

  template< class FunctionType >
  static RangeFieldType norm ( const FunctionType &function,
                               DiscreteFunctionType &discreteFunction,
                               double time )
  {
    const DiscreteFunctionSpaceType &discreteFunctionSpace
      = discreteFunction.space();
    const int polOrder = 2 * discreteFunctionSpace.order() + 2;
    return norm( function, discreteFunction, time, polOrder );
  }
};



// calculates \Vert u - u_h \Vert_{H^1}
template< class DiscreteFunctionImp >
class H1Error
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
    
  typedef typename DiscreteFunctionType :: FunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType :: LocalFunctionType
    LocalFunctionType;
 
  typedef typename DiscreteFunctionSpaceType :: DomainFieldType
    DomainFieldType;
  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType :: RangeFieldType
    RangeFieldType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType :: JacobianRangeType
    JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

  enum { DimDomain = DiscreteFunctionSpaceType :: DimDomain };
  enum { DimRange = DiscreteFunctionSpaceType :: DimRange };

  typedef typename GridPartType :: GridType GridType;
  typedef typename GridType :: template Codim< 0 > :: Entity Entity0Type;
  typedef typename Entity0Type :: Geometry GeometryType;
 
  
public:
  template< class FunctionType >
  static RangeFieldType norm ( const FunctionType &function,
                               DiscreteFunctionType &discreteFunction,
                               double time,
                               int polOrder )
  {
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

    const DiscreteFunctionSpaceType &discreteFunctionSpace
      = discreteFunction.space();

    RangeFieldType error( 0 );
    IteratorType endit = discreteFunctionSpace.end();
    for( IteratorType it = discreteFunctionSpace.begin(); it != endit; ++it )
    {
      CachingQuadrature< GridPartType, 0 > quadrature( *it, polOrder );
      LocalFunctionType localFunction = discreteFunction.localFunction( *it );

      const GeometryType &geometry = (*it).geometry();
      
      const int numQuadraturePoints = quadrature.nop();
      for( int qp = 0; qp < numQuadraturePoints; ++qp )
      {
        const DomainType &x = quadrature.point( qp );
        const DomainType &y = geometry.global( x );
        const double weight
          = quadrature.weight( qp ) * geometry.integrationElement( x );

        RangeType phi, psi;
        function.evaluate( y, time, phi );
        localFunction.evaluate( quadrature[ qp ], psi );

        JacobianRangeType Dphi, Dpsi;
        function.jacobian( y, time, Dphi );
        localFunction.jacobian( quadrature[ qp ], Dpsi );

        for( int i = 0; i < DimRange; ++i ) {
          RangeFieldType localError = SQR( phi[ i ] - psi[ i ] );
          for( int j = 0; j < DimDomain; ++j )
            localError += SQR(Dphi[ i ][ j ] - Dpsi[ i ][ j ]);
          error += weight * localError;
        }
      }
    }
    
    return sqrt( error );
  }

  template< class FunctionType >
  static RangeFieldType norm ( const FunctionType &function,
                               DiscreteFunctionType &discreteFunction,
                               double time )
  {
    const DiscreteFunctionSpaceType &discreteFunctionSpace
      = discreteFunction.space();
    const int polOrder = 2 * discreteFunctionSpace.order() + 2;
    return norm( function, discreteFunction, time, polOrder );
  }
};

 

void adapt ( GridType &grid, DiscreteFunctionType &solution, int step )
{
  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
  
  const DiscreteFunctionSpaceType &discreteFunctionSpace = solution.space();
  
  RestrictProlongOperatorType rp( solution );

  #if GENERIC_ADAPT 
    AdaptationManagerType adaptationManager( grid, rp );
  #else 
    DofManagerType &dofManager = DofManagerFactoryType :: getDofManager( grid );
  #endif  

  std :: string message = (step < 0 ? "Coarsening..." : "Refining..." );
  const int mark = (step < 0 ? -1 : 1);
  const int count = std :: abs( step );
  
  for( int i = 0; i < count; ++i ) {
    IteratorType it = discreteFunctionSpace.begin();
    const IteratorType endit = discreteFunctionSpace.end();
    for( ; it != endit; ++it )
      grid.mark( mark, it );
      
    #if GENERIC_ADAPT
      std :: cout << message << std::endl;
      adaptationManager.adapt();
    #else
      std :: cout << message << " (nongeneric)" << std :: endl;
      grid.adapt( dofManager, rp );
    #endif
  }
}



void algorithm ( GridType &grid, 
                 DiscreteFunctionType &solution, 
                 int step, int turn )
{
  FunctionSpaceType functionSpace;
  ExactSolution f( functionSpace );
  
  LagrangeInterpolation< DiscreteFunctionType > :: interpolateFunction( f, solution );
  double preL2error = L2Error< DiscreteFunctionType > :: norm( f, solution, 0 );
  double preH1error = H1Error< DiscreteFunctionType > :: norm( f, solution, 0 );

  std :: cout << "L2 error before adaption: " << preL2error << std :: endl;
  std :: cout << "H1 error before adaption: " << preH1error << std :: endl; 
  
  adapt( grid, solution, step );
  
  double postL2error = L2Error< DiscreteFunctionType > :: norm( f, solution, 0 );
  double postH1error = H1Error< DiscreteFunctionType > :: norm( f, solution, 0 );

  std :: cout << "L2 error after "
              << (step < 0 ? "restriction" : "prolongation")
              << ": " << postL2error << std :: endl;
  std :: cout << "H1 error after "
              << (step < 0 ? "restriction" : "prolongation")
              << ": " << postH1error << std :: endl; 
  
  #if USE_GRAPE && SHOW_RESTRICT_PROLONG
    if( turn > 0 ) {
      GrapeDataDisplay< GridType > grape( grid );
      grape.dataDisplay( solution );
    }
  #endif

  LagrangeInterpolation< DiscreteFunctionType > :: interpolateFunction( f, solution );
  double newL2error = L2Error< DiscreteFunctionType > :: norm( f, solution, 0 );
  double newH1error = H1Error< DiscreteFunctionType > :: norm( f, solution, 0 );

  std :: cout << "L2 error for interpolation after adaption: " << newL2error << std :: endl;
  std :: cout << "H1 error for interpolation after adaption: " << newH1error << std :: endl; 
  
  #if USE_GRAPE && SHOW_INTERPOLATION
    if( turn > 0 ) {
      GrapeDataDisplay< GridType > grape( grid );
      grape.dataDisplay( solution );
    }
  #endif

  double l2eoc = -log( newL2error / preL2error) / M_LN2;
  double h1eoc = -log( newH1error / preH1error) / M_LN2;

  std :: cout << "L2 EOC: " << l2eoc << std :: endl;
  std :: cout << "H1 EOC: " << h1eoc << std :: endl;

  #if WRITE_DATA
    GrapeDataIO< GridType > dataio; 
    dataio.writeGrid( grid, xdr, "gridout", 0, turn );
    dataio.writeData( solution, xdr, "sol", turn );
  #endif

  std :: cout << std :: endl;
}



int main ( int argc, char **argv )
{
  if( argc != 2 ) {
    std :: cerr << "Usage: " << argv[ 0 ] << "<maxlevel>" << std :: endl;
    exit( 1 );
  }
  
  int ml = atoi( argv[ 1 ] );

  char tmp[ 100 ]; 
  sprintf( tmp, "%ddgrid.dgf", GRIDDIM );
  GridPtr< GridType > gridptr( tmp );

  const int step = DGFGridInfo< GridType > :: refineStepsForHalf();

  GridPartType gridPart( *gridptr );
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
  DiscreteFunctionType solution( "solution", discreteFunctionSpace );
  solution.clear();

  std :: cout << std :: endl << "Refining: " << std :: endl;
  for( int i = 0; i < ml; ++i )
    algorithm( *gridptr, solution, step, (i == ml-1) );
  
  std :: cout << std :: endl << "Coarsening:" << std::endl;
  for( int i = ml - 1; i >= 0; --i )
    algorithm( *gridptr, solution, -step, 1 );

  return 0;
}

