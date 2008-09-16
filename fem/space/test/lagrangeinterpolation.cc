#include <config.h>

// Preprocessor Definitions
// ------------------------

// allow higher order Lagrange mappers
#define USE_TWISTFREE_MAPPER

// unset to use L^2 error
#define USE_H1ERROR



// Includes
// --------

#include <iostream>

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/quadrature/cachequad.hh>

#include <dune/fem/operator/lagrangeinterpolation.hh>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh> 

#include <dune/grid/common/referenceelements.hh>

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

using namespace Dune;

// polynom approximation order of quadratures, 
// at least poolynom order of basis functions 
#ifdef POLORDER
  const int polOrder = POLORDER;
#else
  const int polOrder = 1;
#endif

//***********************************************************************
/*! L2 Projection of a function f: 

  This is an example how to solve the equation on 
  \f[\Omega = (0,1)^2 \f]

  \f[ \int_{\Omega} u \phi = \int_{\Omega} f \phi  \ \ \ in \Omega \f]
  \f[ f(x,y) = x ( 1 - x) y ( 1 - y ) \f]

  Here u is the L_2 projection of f. 

  The Projection should converge to the given function f.
  with the finite element method using lagrangian elements of polynom order +1.
*/
//***********************************************************************

//! the index set we are using 
//typedef HierarchicGridPart< GridType > GridPartType;
typedef LeafGridPart< GridType > GridPartType;

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
//typedef MatrixFunctionSpace< double, double, dimworld, 3, 5 >
//  FunctionSpaceType;
typedef FunctionSpace< double, double, dimworld, 1 > FunctionSpaceType;
typedef MatrixFunctionSpace< double, double, dimworld, 1, dimworld >
  GradientFunctionSpaceType;

//! define the function space our unkown belong to 
typedef DiscontinuousGalerkinSpace< GradientFunctionSpaceType, 
                                    GridPartType,
                                    polOrder,
                                    CachingStorage >
  DiscreteGradientFunctionSpaceType;
typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
  DiscreteFunctionSpaceType;

//! define the type of discrete function we are using , see
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
  DiscreteFunctionType;
typedef AdaptiveDiscreteFunction< DiscreteGradientFunctionSpaceType >
  DiscreteGradientFunctionType;

//! Get the Dofmanager type
typedef DofManager< GridType > DofManagerType;
typedef DofManagerFactory< DofManagerType > DofManagerFactoryType;



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
                      
 
// ********************************************************************
template <class DiscreteFunctionType>
class L2Projection
{
  typedef typename DiscreteFunctionType::FunctionSpaceType DiscreteFunctionSpaceType;

 public:
  template <class FunctionType>
  static void project (const FunctionType &f, DiscreteFunctionType &discFunc, int polOrd) 
  {
    typedef typename DiscreteFunctionSpaceType::Traits::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::Traits::IteratorType Iterator;

    const DiscreteFunctionSpaceType& space =  discFunc.space();

    discFunc.clear();

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

    typename DiscreteFunctionSpaceType::RangeType ret (0.0);
    typename DiscreteFunctionSpaceType::RangeType phi (0.0);

    Iterator endit = space.end();
    for(Iterator it = space.begin(); it != endit ; ++it) 
    {
      // Get quadrature rule
      CachingQuadrature<GridPartType,0> quad(*it, polOrd);

      LocalFuncType lf = discFunc.localFunction(*it);

      //! Note: BaseFunctions must be ortho-normal!!!!
      typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType ; 
      const BaseFunctionSetType & baseset =
        lf.baseFunctionSet();

      //const typename GridType::template Codim<0>::Entity::Geometry& 
      //  itGeom = (*it).geometry();
     
      const int quadNop = quad.nop();
      const int numDofs = lf.numDofs();
      for(int qP = 0; qP < quadNop ; ++qP) 
      {
        // f.evaluate(itGeom.global(quad.point(qP)), ret);
        f.localFunction(*it).jacobian(quad,qP, ret);
        for(int i=0; i<numDofs; ++i) {
          baseset.evaluate(i,quad[qP],phi);
          lf[i] += quad.weight(qP) * (ret * phi) ;
        }
      }
    }
  }
  
  template <class FunctionType>
  static void project (const FunctionType &f, DiscreteFunctionType &discFunc) 
  {
    const DiscreteFunctionSpaceType& space =  discFunc.space();
    int polOrd = 2 * space.order();
    project(f,discFunc,polOrd);
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
  
  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
  
  enum { DimRange = DiscreteFunctionSpaceType :: DimRange };

  typedef typename GridPartType :: GridType GridType;
  typedef typename GridType :: template Codim< 0 > :: Entity Entity0Type;
  typedef typename Entity0Type :: Geometry GeometryType;
 
  
public:
  template< class FunctionType >
  static RangeType norm ( const FunctionType &function,
                          DiscreteFunctionType &discreteFunction,
                          double time,
                          int polOrder )
  {
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

    const DiscreteFunctionSpaceType &discreteFunctionSpace
      = discreteFunction.space();

    RangeType error( 0 );
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
        localFunction.evaluate( quadrature[qp], psi );

        for( int i = 0; i < DimRange; ++i )
          error[ i ] += weight * SQR( phi[ i ] - psi[ i ] );
      }
    }
    
    for( int i = 0; i < DimRange; ++i ) 
      error[ i ] = sqrt( error[ i ] );
    
    return error;
  }

  template< class FunctionType >
  static RangeType norm ( const FunctionType &function,
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
  static RangeType norm ( const FunctionType &function,
                          DiscreteFunctionType &discreteFunction,
                          double time,
                          int polOrder )
  {
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;

    const DiscreteFunctionSpaceType &discreteFunctionSpace
      = discreteFunction.space();

    RangeType error( 0 );
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
        localFunction.evaluate( quadrature[qp], psi );

        JacobianRangeType Dphi, Dpsi;
        function.jacobian( y, time, Dphi );
        localFunction.jacobian( quadrature[ qp ], Dpsi );

        for( int i = 0; i < DimRange; ++i ) {
          RangeFieldType localError = SQR( phi[ i ] - psi[ i ] );
          for( int j = 0; j < DimDomain; ++j )
            localError += SQR(Dphi[ i ][ j ] - Dpsi[ i ][ j ]);
          error[ i ] += weight * localError;
        }
      }
    }
    
    for( int i = 0; i < DimRange; ++i ) 
      error[ i ] = sqrt( error[ i ] );
    
    return error;
  }

  template< class FunctionType >
  static RangeType norm ( const FunctionType &function,
                          DiscreteFunctionType &discreteFunction,
                          double time )
  {
    const DiscreteFunctionSpaceType &discreteFunctionSpace
      = discreteFunction.space();
    const int polOrder = 2 * discreteFunctionSpace.order() + 2;
    return norm( function, discreteFunction, time, polOrder );
  }
};



double algorithm( GridType &grid, DiscreteFunctionType &solution, int turn )
{
  typedef DiscreteFunctionSpaceType :: RangeType RangeType; 

  GridPartType part( grid );
  DiscreteFunctionSpaceType discreteFunctionSpace( part );
  ExactSolution f( discreteFunctionSpace );
  
  //! perform Lagrange interpolation
  LagrangeInterpolation< DiscreteFunctionType >
  :: interpolateFunction( f, solution );

  #if 0
  DiscreteGradientFunctionSpaceType discreteGradientFunctionSpace( part );
  DiscreteGradientFunctionType graddf("grad", discreteGradientFunctionSpace );
  
  //! perform l2-projection
  L2Projection< DiscreteGradientFunctionType > :: project( solution, graddf );
  #endif

  // calculation L2 error 
  // pol ord for calculation the error chould by higher than 
  // pol for evaluation the basefunctions
  #ifdef USE_H1ERROR
    RangeType error = H1Error< DiscreteFunctionType > :: norm( f ,solution, 0 );
    std :: cout << std :: endl << "H1 Error: [ " << error[ 0 ];
  #else
    RangeType error = L2Error< DiscreteFunctionType > :: norm( f ,solution, 0 );
    std :: cout << std :: endl << "L2 Error: [ " << error[ 0 ];
  #endif
  for( int i = 1; i < RangeType :: dimension; ++i )
    std :: cout << ", " << error[ i ];
  std :: cout << " ]" << std :: endl;

  #if HAVE_GRAPE
   // if Grape was found, then display last solution 
   if(turn > 0)
   {
     GrapeDataDisplay < GridType > grape(part); 
     grape.addData( solution );
     //grape.addData( graddf );
     grape.display( );
   }
  #endif
   
  return sqrt( error * error );
}


//**************************************************
//
//  main programm, run algorithm twice to calc EOC 
//
//**************************************************
int main (int argc, char **argv)
{
  std :: cout << "L2 projection test for polynomial order " << polOrder
              << std :: endl;
    
  if( argc != 2)
  {
    std :: cerr << "Usage: " << argv[ 0 ] << " <maxlevel>" << std :: endl;
    exit( 1 );
  }

  try
  {
    int ml = atoi( argv[1] );
    double error[ ml ];

    char tmp[16];
    sprintf( tmp, "%d", dimworld );
    std :: string macroGridName( tmp );
    macroGridName += "dgrid.dgf";

    GridPtr<GridType> gridptr(macroGridName);
    GridType& grid=*gridptr;
    const int step = Dune::DGFGridInfo<GridType>::refineStepsForHalf();
    GridPartType gridPart( grid );

    for( int i = 0; i < ml; ++i )
    {
      grid.globalRefine( step );
      DofManagerType& dm = DofManagerFactoryType :: getDofManager( grid );
      dm.resize();

      DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
      DiscreteFunctionType solution( "sol", discreteFunctionSpace );
      solution.clear();
      
      error[ i ] = algorithm( grid, solution, i == ml-1 );
      if( i > 0 ) {
        double eoc = log( error[ i-1 ] / error[ i ] ) / M_LN2;
        std :: cout << "EOC = " << eoc << std :: endl;
      }
    }
    return 0;
  }
  catch( Exception e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}

