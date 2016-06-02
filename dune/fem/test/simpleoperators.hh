#ifndef MASSOPERATOR_HH
#define MASSOPERATOR_HH

#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>

template< class DomainFunction, class RangeFunction = DomainFunction >
struct NullOp
  : public Dune::Fem::Operator< DomainFunction, RangeFunction >
{
  template< class ... Args >
  NullOp ( Args ... args ) {}
  NullOp ( std::string &name, int value ) : name_( name )
  {
    std::cout<<value <<std::endl;
  }
  void operator() ( const DomainFunction &arg, RangeFunction &dest ) const { dest.clear(); }

  std::string name_;
};


template< class DomainFunction, class RangeFunction = DomainFunction >
class SimpleMassOperator
  : public Dune::Fem::Operator< DomainFunction, RangeFunction >
{
  typedef SimpleMassOperator< DomainFunction, RangeFunction > ThisType;
  typedef Dune::Fem::Operator< DomainFunction, RangeFunction > BaseType;

public:
  typedef typename BaseType::DomainFunctionType DomainFunctionType;
  typedef typename BaseType::RangeFunctionType RangeFunctionType;

  typedef typename DomainFunctionType::GridPartType GridPartType;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  std::string name_;

  template< class ... Args >
  SimpleMassOperator ( Args ... args ) {}

  SimpleMassOperator ( std::string &name, int value ) : name_( name )
  {
    std::cout<<value <<std::endl;
  }

  void operator() ( const DomainFunctionType &arg, RangeFunctionType &dest ) const
  {
    dest.clear();
    Dune::Fem::ConstLocalFunction< DomainFunctionType > argLocal( arg );
    Dune::Fem::TemporaryLocalFunction< typename RangeFunctionType::DiscreteFunctionSpaceType > local( dest.space() );

    for( const auto &entity : arg.space() )
    {
      auto geometry = entity.geometry();

      argLocal.init( entity );
      local.init( entity );

      // run over quadrature points
      QuadratureType quadrature( entity, 2* argLocal.order()+1 );
      for( const auto qp : quadrature )
      {
        // evaluate u
        const typename QuadratureType::CoordinateType &x = qp.position();

        typename DomainFunctionType::RangeType uValue;
        argLocal.evaluate( qp, uValue );

        // put all things together and don't forget quadrature weights
        const double weight = qp.weight()*geometry.integrationElement( x );

        // apply weight
        uValue *= weight;

        // add to local function
        local.axpy( qp, uValue );
      }
      dest.addLocalDofs( entity, local.localDofVector() );
    }
    dest.communicate();
  }
};



template< class DomainFunction, class RangeFunction = DomainFunction >
class SimpleLaplaceOperator
  : public Dune::Fem::Operator< DomainFunction, RangeFunction >
{
  typedef SimpleLaplaceOperator< DomainFunction, RangeFunction > ThisType;
  typedef Dune::Fem::Operator< DomainFunction, RangeFunction > BaseType;

public:
  typedef typename BaseType::DomainFunctionType DomainFunctionType;
  typedef typename BaseType::RangeFunctionType RangeFunctionType;

  typedef typename DomainFunctionType::GridPartType GridPartType;

  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  std::string name_;

  template< class ... Args >
  SimpleLaplaceOperator ( Args ... args ) {}

  template< class ... Args >
  SimpleLaplaceOperator ( std::string &name, Args ... args ) : name_( name ) {}

  void operator() ( const DomainFunctionType &arg, RangeFunctionType &dest ) const
  {
    dest.clear();
    Dune::Fem::ConstLocalFunction< DomainFunctionType > argLocal( arg );
    Dune::Fem::TemporaryLocalFunction< typename RangeFunctionType::DiscreteFunctionSpaceType > local( dest.space() );

    for( const auto &entity : arg.space() )
    {
      auto geometry = entity.geometry();

      argLocal.init( entity );
      local.init( entity );

      // run over quadrature points
      QuadratureType quadrature( entity, 2* argLocal.order()+1 );
      for( const auto qp : quadrature )
      {
        // evaluate u
        const typename QuadratureType::CoordinateType &x = qp.position();

        typename DomainFunctionType::JacobianRangeType uValue;
        argLocal.jacobian( qp, uValue );

        // put all things together and don't forget quadrature weights
        const double weight = qp.weight()*geometry.integrationElement( x );

        // apply weight
        uValue *= weight;

        // add to local function
        local.axpy( qp, uValue );
      }
      dest.addLocalDofs( entity, local.localDofVector() );
    }
    dest.communicate();
  }
};

#endif // #ifndef MASSOPERATOR_HH
