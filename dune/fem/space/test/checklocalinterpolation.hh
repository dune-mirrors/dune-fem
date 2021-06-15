#ifndef DUNE_FEM_TEST_CHECKLOCALINTERPOLATION_HH
#define DUNE_FEM_TEST_CHECKLOCALINTERPOLATION_HH

#include <iostream>

#include <dune/common/classname.hh>
#include <dune/common/dynvector.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/common/localinterpolation.hh>

template< class DiscreteSpace >
struct LocalBasis
{
  typedef DiscreteSpace FunctionSpaceType;
  typedef typename DiscreteSpace::BasisFunctionSetType BasisFunctionSetType;
  typedef typename DiscreteSpace::DomainType DomainType;
  typedef typename DiscreteSpace::RangeType RangeType;
  typedef typename DiscreteSpace::JacobianRangeType JacobianRangeType;
  typedef typename DiscreteSpace::HessianRangeType HessianRangeType;
  typedef typename DiscreteSpace::RangeFieldType RangeFieldType;

  typedef typename BasisFunctionSetType::EntityType EntityType;

  LocalBasis ( const BasisFunctionSetType &bSet, int i )
    : bSet_( bSet ),
    comp_( bSet.size(), 0 )
  {
    comp_[ i ] = 1.0;
  }

  template< class Point >
  void evaluate ( const Point &arg, RangeType &dest ) const
  {
    dest = 0;
    bSet_.evaluateAll( arg, comp_, dest );
  }
  template <class Point>
  RangeType operator()(const Point &arg) const
  {
    RangeType dest;
    evaluate(arg,dest);
    return dest;
  }


  template< class Point >
  RangeType evaluate ( const Point &arg ) const
  {
    RangeType ret;
    evaluate( arg, ret );
    return ret;
  }

  template< class Point >
  void jacobian ( const Point &arg, JacobianRangeType &jac ) const
  {
    jac = 0;
    bSet_.jacobianAll( arg, comp_, jac );
  }

  template< class Point >
  JacobianRangeType jacobian ( const Point &arg ) const
  {
    JacobianRangeType jac;
    jacobian( arg, jac );
    return jac;
  }

  template< class QuadratureType, class Vector >
  void evaluateQuadrature( const QuadratureType &quad, Vector& vec ) const
  {
    evaluateQuadrature( quad, vec, vec[ 0 ] );
  }

  int order () const { return bSet_.order(); }
  const EntityType &entity () const { return bSet_.entity(); }

private:
  template< class QuadratureType, class Vector >
  void evaluateQuadrature( const QuadratureType &quad, Vector& vec, const RangeType&  ) const
  {
    const int nop = quad.nop();
    for( int qp = 0; qp < nop; ++qp )
    {
      evaluate( quad[ qp ], vec[ qp ] );
    }
  }

  template< class QuadratureType, class Vector >
  void evaluateQuadrature( const QuadratureType &quad, Vector& vec, const JacobianRangeType&  ) const
  {
    const int nop = quad.nop();
    for( int qp = 0; qp < nop; ++qp )
    {
      jacobian( quad[ qp ], vec[ qp ] );
    }
  }

  const BasisFunctionSetType &bSet_;
  Dune::DynamicVector< RangeFieldType > comp_;
};


template< class Space >
auto checkLocalInterpolation ( const Space &space )
  -> std::enable_if_t< Dune::Fem::Capabilities::hasInterpolation< Space >::v >
{
  std::cout << ">>> Testing local interpolation of " << Dune::className< Space >() << ":" << std::endl;
  typedef LocalBasis< Space > LocalBasisType;

  Dune::Fem::LocalInterpolation< Space > interpolation( space );

  for( auto entity : space )
  {
    auto bSet = space.basisFunctionSet( entity );
    auto guard = Dune::Fem::bindGuard( interpolation, entity );

    for( std::size_t i = 0; i < bSet.size(); ++i )
    {
      LocalBasisType local( bSet, i );

      std::vector< typename Space::RangeFieldType > phii( bSet.size(), 0 );
      interpolation( local, phii );

      for( std::size_t j = 0; j < phii.size(); ++j )
        if( std::abs( phii[ j ] -( i == j ))  > 1e-8 )
        {
          std::cerr << "Local interpolation error for basis function: "<<j <<" with diff: " <<  phii[ j ] <<  std::endl;
          // std::abort();
        }
    }
  }
  std::cout <<"... done." <<std::endl;
}

template< class Space >
auto checkLocalInterpolation ( const Space &space )
  -> std::enable_if_t< !Dune::Fem::Capabilities::hasInterpolation< Space >::v >
{
  std::cout << "Nothing to be done for spaces without local interpolation" <<std::endl;
}

#endif // #ifndef DUNE_FEM_TEST_CHECKLOCALINTERPOLATION_HH
