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

// #include <dgfgridtype.hh>

#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/space/common/interpolate.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include <dune/fem/test/testgrid.hh>

using namespace Dune;
using namespace Fem;

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
typedef GridSelector::GridType MyGridType;

//! the index set we are using
//typedef HierarchicGridPart< MyGridType > GridPartType;
typedef LeafGridPart< MyGridType > GridPartType;

//! define the function space our unkown belong to
typedef DiscontinuousGalerkinSpace< Dune::FieldMatrix< double, 1, MyGridType::dimensionworld >, GridPartType, polOrder > DiscreteGradientFunctionSpaceType;
typedef LagrangeDiscreteFunctionSpace< Dim< 1 >, GridPartType, polOrder > DiscreteFunctionSpaceType;

//! define the type of discrete function we are using , see
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
  DiscreteFunctionType;
typedef AdaptiveDiscreteFunction< DiscreteGradientFunctionSpaceType >
  DiscreteGradientFunctionType;

//! Get the Dofmanager type
typedef DofManager< MyGridType > DofManagerType;

class ExactSolution
: public Fem::Function< typename DiscreteFunctionSpaceType::FunctionSpaceType, ExactSolution >
{
  typedef ExactSolution ThisType;
  typedef Fem::Function< FunctionSpaceType, ExactSolution > BaseType;

public:
  typedef BaseType::DomainFieldType DomainFieldType;
  typedef BaseType::RangeFieldType RangeFieldType;

  typedef BaseType::DomainType DomainType;
  typedef BaseType::RangeType RangeType;
  typedef BaseType::JacobianRangeType JacobianRangeType;

  void evaluate ( const DomainType &x, RangeType &phi ) const
  {
    phi = 1;
    for( int i = 0; i < DomainType :: dimension; ++i )
      // phi[ 0 ] += x[ i ] * x[ i ];
      phi[ 0 ] *= sin( M_PI * x[ i ] );
  }

  void jacobian( const DomainType &x, JacobianRangeType &Dphi ) const
  {
    Dphi = 1;
    for( int i = 0; i < DomainType :: dimension; ++i )
      for( int j = 0; j < DomainType :: dimension; ++j )
        // Dphi[ 0 ][ j ] *= ((i != j) ? 1. : 2.*x[i]);
        Dphi[ 0 ][ j ] *= ((i != j) ? sin( M_PI * x[ i ]) : M_PI * cos( M_PI * x[ i ] ));
  }

  void evaluate( const DomainType &x, RangeFieldType t, RangeType &phi ) const
  {
    evaluate( x, phi );
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

  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

  enum { DimRange = DiscreteFunctionSpaceType :: DimRange };

  typedef typename GridPartType::GridType GridType;
  typedef typename GridType::template Codim< 0 >::Entity Entity0Type;
  typedef typename Entity0Type::Geometry GeometryType;


public:
  template< class FunctionType >
  static RangeType norm ( const FunctionType &function,
                          DiscreteFunctionType &discreteFunction,
                          double time,
                          int polOrder )
  {
    const DiscreteFunctionSpaceType &discreteFunctionSpace
      = discreteFunction.space();

    Dune::Fem::ConstLocalFunction< DiscreteFunctionType > dfLocal( discreteFunction );

    RangeType error( 0 );
    for( const auto& entity : discreteFunctionSpace )
    {
      CachingQuadrature< GridPartType, 0 > quadrature( entity, polOrder );
      auto guard = Dune::Fem::bindGuard( dfLocal, entity );

      const GeometryType &geometry = entity.geometry();

      const int numQuadraturePoints = quadrature.nop();
      for( int qp = 0; qp < numQuadraturePoints; ++qp )
      {
        const DomainType &x = quadrature.point( qp );
        const double weight
          = quadrature.weight( qp ) * geometry.integrationElement( x );

        RangeType phi, psi;
        function.evaluate( geometry.global( x ), time, phi );
        dfLocal.evaluate( quadrature[qp], psi );

        for( int i = 0; i < DimRange; ++i )
          error[ i ] += weight * ((phi[ i ] - psi[ i ])*(phi[ i ] - psi[ i ]));
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

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  static const int DimDomain = DiscreteFunctionSpaceType::DimDomain;
  static const int DimRange = DiscreteFunctionSpaceType::DimRange;

  typedef typename GridPartType::GridType GridType;
  typedef typename GridType::template Codim< 0 >::Entity Entity0Type;
  typedef typename Entity0Type::Geometry GeometryType;


public:
  template< class FunctionType >
  static RangeType norm ( const FunctionType &function,
                          DiscreteFunctionType &discreteFunction,
                          double time,
                          int polOrder )
  {
    const DiscreteFunctionSpaceType &discreteFunctionSpace
      = discreteFunction.space();

    Dune::Fem::ConstLocalFunction< DiscreteFunctionType > dfLocal( discreteFunction );

    RangeType error( 0 );
    for( const auto& entity : discreteFunctionSpace )
    {
      CachingQuadrature< GridPartType, 0 > quadrature( entity, polOrder );
      auto guard = Dune::Fem::bindGuard( dfLocal, entity );

      const GeometryType &geometry = entity.geometry();

      const int numQuadraturePoints = quadrature.nop();
      for( int qp = 0; qp < numQuadraturePoints; ++qp )
      {
        const DomainType &x = quadrature.point( qp );
        const DomainType &y = geometry.global( x );
        const double weight
          = quadrature.weight( qp ) * geometry.integrationElement( x );

        RangeType phi, psi;
        function.evaluate( y, time, phi );
        dfLocal.evaluate( quadrature[qp], psi );

        JacobianRangeType Dphi, Dpsi;
        function.jacobian( y, time, Dphi );
        dfLocal.jacobian( quadrature[ qp ], Dpsi );

        for( int i = 0; i < DimRange; ++i ) {
          RangeFieldType localError = (phi[ i ] - psi[ i ])*(phi[ i ] - psi[ i ]);
          for( int j = 0; j < DimDomain; ++j )
            localError += (Dphi[ i ][ j ] - Dpsi[ i ][ j ])*(Dphi[ i ][ j ] - Dpsi[ i ][ j ]);
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



double algorithm( MyGridType &grid, DiscreteFunctionType &solution, int turn )
{
  typedef DiscreteFunctionSpaceType :: RangeType RangeType;

  GridPartType part( grid );
  DiscreteFunctionSpaceType discreteFunctionSpace( part );
  ExactSolution f;

  //! perform Lagrange interpolation
  interpolate( gridFunctionAdapter( f, part, discreteFunctionSpace.order() + 2 ), solution );

  // calculation L2 error
  // pol ord for calculation the error chould by higher than
  // pol for evaluation the basefunctions
  #ifdef USE_H1ERROR
    RangeType error = H1Error< DiscreteFunctionType >::norm( f ,solution, 0 );
    std :: cout << std :: endl << "H1 Error: [ " << error[ 0 ];
  #else
    RangeType error = L2Error< DiscreteFunctionType >::norm( f ,solution, 0 );
    std :: cout << std :: endl << "L2 Error: [ " << error[ 0 ];
  #endif
  for( int i = 1; i < RangeType :: dimension; ++i )
    std :: cout << ", " << error[ i ];
  std :: cout << " ]" << std :: endl;

  #if HAVE_GRAPE
   // if Grape was found, then display last solution
   if(turn > 0)
   {
     GrapeDataDisplay< MyGridType > grape(part);
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
try
{
  MPIManager::initialize( argc, argv );

  std :: cout << "L2 projection test for polynomial order " << polOrder
              << std :: endl;

  int ml = 2;
  if( argc > 1 )
    ml = atoi( argv[1] );

  std::vector< double > error( ml );

  MyGridType &grid = Dune::Fem::TestGrid::grid();
  const int step = Dune::Fem::TestGrid::refineStepsForHalf();
  GridPartType gridPart( grid );

  for( int i = 0; i < ml; ++i )
  {
    Dune::Fem::GlobalRefine::apply( grid, step );

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
catch( const Dune::Exception &e )
{
  std::cerr << "Error: " << e << std::endl;
  return 1;
}
