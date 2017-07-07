#include <config.h>

#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

static const int dimw = Dune::GridSelector::dimworld;

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/lagrange.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>

#include <dune/fem/misc/double.hh>

#include <dune/fem/test/testgrid.hh>


// polynom approximation order of quadratures,
// at least polynom order of basis functions
const int polOrd = POLORDER;

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
typedef Dune::GridSelector::GridType MyGridType;
typedef Dune::Fem::AdaptiveLeafGridPart< MyGridType > GridPartType;

// see dune/common/functionspace.hh
typedef Dune::Fem::FunctionSpace< MyGridType::ctype, double, dimw, 2 > FuncSpace1;
typedef Dune::Fem::FunctionSpace< MyGridType::ctype, double, dimw, 1 > FuncSpace2;

typedef FuncSpace1::RangeType RangeType1;
typedef FuncSpace2::RangeType RangeType2;
typedef FuncSpace1::JacobianRangeType JacobianRangeType1;
typedef FuncSpace2::JacobianRangeType JacobianRangeType2;
typedef FuncSpace1::HessianRangeType HessianRangeType1;
typedef FuncSpace2::HessianRangeType HessianRangeType2;

//! define the function space our unkown belong to
//! see dune/fem/lagrangebase.hh
typedef Dune::Fem::LagrangeDiscontinuousGalerkinSpace< FuncSpace1, GridPartType,
                                            polOrd+1, Dune::Fem::CachingStorage > DiscreteFunctionSpaceType1;

typedef Dune::Fem::LagrangeDiscontinuousGalerkinSpace< FuncSpace2, GridPartType,
                                            polOrd, Dune::Fem::CachingStorage > DiscreteFunctionSpaceType2;

typedef Dune::Fem::TupleDiscreteFunctionSpace< DiscreteFunctionSpaceType1, DiscreteFunctionSpaceType2 > DiscreteFunctionSpaceType;

typedef DiscreteFunctionSpaceType::IteratorType::Entity EntityType;

typedef DiscreteFunctionSpaceType::FunctionSpaceType FuncSpace;
typedef FuncSpace::RangeType RangeType;
typedef FuncSpace::JacobianRangeType JacobianRangeType;
typedef FuncSpace::HessianRangeType HessianRangeType;


const int dimRange = FuncSpace::dimRange;
const int dimRange1 = FuncSpace1::dimRange;
const int dimRange2 = FuncSpace2::dimRange;
const int dimDomain = FuncSpace::dimDomain;


void checkBasisSet ( const EntityType &entity, const DiscreteFunctionSpaceType &space )
{
  auto bSet = space.basisFunctionSet( entity );

  auto bSet1 = space.subDiscreteFunctionSpace< 0 >().basisFunctionSet( entity );
  auto bSet2 = space.subDiscreteFunctionSpace< 1 >().basisFunctionSet( entity );

  const int nrBasis = bSet.size();
  const int nrBasis1 = bSet1.size();
  const int nrBasis2 = bSet2.size();

  std::vector< RangeType > ranges( nrBasis, RangeType( 0. ) );
  std::vector< RangeType1 > ranges1( nrBasis1, RangeType1( 0. ) );
  std::vector< RangeType2 > ranges2( nrBasis2, RangeType2( 0. ) );

  std::vector< JacobianRangeType > jacs( nrBasis, JacobianRangeType( 0. ) );
  std::vector< JacobianRangeType1 > jacs1( nrBasis1, JacobianRangeType1( 0. ) );
  std::vector< JacobianRangeType2 > jacs2( nrBasis2, JacobianRangeType2( 0. ) );

  HessianRangeType null;
  HessianRangeType1 null1;
  HessianRangeType2 null2;

  for( int i = 0; i < dimDomain; ++i )
    for( int j = 0; j < dimDomain; ++j )
    {
      for( int r = 0; r < dimRange; ++r )
        null[ r ][ i ][ j ] = .0;
      for( int r = 0; r < dimRange1; ++r )
        null1[ r ][ i ][ j ] = .0;
      for( int r = 0; r < dimRange2; ++r )
        null2[ r ][ i ][ j ] = .0;
    }

  std::vector< HessianRangeType > hessians( nrBasis, null );
  std::vector< HessianRangeType1 > hessians1( nrBasis1, null1 );
  std::vector< HessianRangeType2 > hessians2( nrBasis2, null2 );

  Dune::Fem::CachingQuadrature< GridPartType, 0 > quadrature( entity, space.order() + 1 );
  for( const auto& qp : quadrature )
  {

    bSet.evaluateAll( qp, ranges );
    bSet1.evaluateAll( qp, ranges1 );
    bSet2.evaluateAll( qp, ranges2 );

    bSet.jacobianAll( qp, jacs );
    bSet1.jacobianAll( qp, jacs1 );
    bSet2.jacobianAll( qp, jacs2 );

    bSet.hessianAll( qp, hessians );
    bSet1.hessianAll( qp, hessians1 );
    bSet2.hessianAll( qp, hessians2 );

    RangeType value( 0 );
    JacobianRangeType jac( 0 );
    HessianRangeType hessian( null );

    for( int i = 0; i < nrBasis; ++i )
    {
      for( int r = 0; r < dimRange; ++r )
      {
        if( r < dimRange1 )
        {
          for( int j = 0; j < dimDomain; ++j )
          {
            jac[ r ][ j ] = ( i < nrBasis1 ) ? jacs1[ i ][ r ][ j ] : 0.;

            for( int k = 0; k < dimDomain; ++k )
              hessian[ r ][ j ][ k ] = ( i < nrBasis1 ) ? hessians1[ i ][ r ][ j ][ k ] : 0.;
          }

          value[ r ] = ( i < nrBasis1 ) ? ranges1[ i ][ r ] : 0.;
        }
        else
        {
          value[ r ] = ( i < nrBasis1 ) ? 0. : ranges2[ i - nrBasis1 ][ r - dimRange1 ];

          for( int j = 0; j < dimDomain; ++j )
          {
            jac[ r ][ j ] = ( i < nrBasis1 ) ? 0. : jacs2[ i-nrBasis1 ][ r - dimRange1 ][ j ];

            for( int k = 0; k < dimDomain; ++k )
              hessian[ r ][ j ][ k ] = ( i < nrBasis1 ) ? 0. : hessians2[ i - nrBasis1 ][ r - dimRange1 ][ j ][ k ];
          }
        }
      }

      // check function value
      value -= ranges[ i ];
      if( value.two_norm() > 1e-8 )
        DUNE_THROW( Dune::InvalidStateException, "Basisfunction::evaluate returns wrong value." );

      JacobianRangeType jacHelp( jac );
      // check jac
      jacHelp -= jacs[ i ];
      if( jacHelp.frobenius_norm() > 1e-8 )
        DUNE_THROW( Dune::InvalidStateException, "Basisfunction::jacobian returns wrong value." );

      for( int r = 0; r < dimRange; ++r )
        hessian[ r ] -= hessians[ i ][ r ];

      double val = 0;
      for( int r = 0; r < dimRange; ++r )
        val += hessian[ r ].frobenius_norm2();

      if( std::sqrt( val ) > 1e-8 )
        DUNE_THROW( Dune::InvalidStateException, "Basisfunction::hessian returns wrong value." );
    }
  }
}


struct DummyLocalFunction
{
  typedef FuncSpace FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  static const int dimRange = FunctionSpaceType::dimRange;
  static const int dimDomain = FunctionSpaceType::dimDomain;

  typedef typename DiscreteFunctionSpaceType::EntityType EntityType;

  DummyLocalFunction ( const EntityType &entity ) : entity_( entity ){}

  template< class Point >
  void evaluate ( const Point &arg, RangeType &dest ) const
  {
    dest = 1.;
  }

  template< class Point >
  void jacobian ( const Point &arg, JacobianRangeType &jac ) const
  {
    jac = 0;
  }

  template< class Point >
  void hessian ( const Point &arg, HessianRangeType &hes ) const
  {}

  int order () const { return 1; }
  const EntityType &entity () const { return entity_; }

  void init ( const EntityType &entity ) {}

protected:
  const EntityType &entity_;
};


void checkInterpolation ( const EntityType &entity, const DiscreteFunctionSpaceType &space )
{
  std::vector< double > ldv;
  ldv.resize( space.basisFunctionSet( entity ).size() );

  DummyLocalFunction lf( entity );
  space.interpolation( entity ) ( lf, ldv );

  //TODO check wether interpolation is successfull
}


int main ( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );
  try
  {
    MyGridType &grid = Dune::Fem::TestGrid::grid();

    GridPartType part( grid );
    DiscreteFunctionSpaceType space( part );

    EntityType entity = *(++space.begin());

    checkBasisSet( entity, space );
    checkInterpolation( entity, space );

    return 0;
  }
  catch( const Dune::Exception &exception )
  {
    std::cerr << exception << std::endl;
    return 1;
  }
}

