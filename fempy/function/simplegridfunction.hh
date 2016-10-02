#ifndef DUNE_FEMPY_FUNCTION_SIMPLEGRIDFUNCTION_HH
#define DUNE_FEMPY_FUNCTION_SIMPLEGRIDFUNCTION_HH

#include <cassert>

#include <limits>
#include <type_traits>
#include <utility>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/space/common/functionspace.hh>

namespace Dune
{

  namespace FemPy
  {

    // SimpleLocalFunction
    // -------------------

    template< class GridPart, class LocalEvaluator >
    class SimpleLocalFunction
    {
      typedef SimpleLocalFunction< GridPart, LocalEvaluator > This;

    public:
      typedef typename GridPart::template Codim< 0 >::EntityType EntityType;

      typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;
      typedef std::decay_t< std::result_of_t< LocalEvaluator( EntityType, LocalCoordinateType ) > > Value;

      typedef Fem::FunctionSpace< typename GridPart::ctype, typename FieldTraits< Value >::field_type, GridPart::dimensionworld, Value::dimension > FunctionSpaceType;

      static const int dimDomain = FunctionSpaceType::dimDomain;
      static const int dimRange = FunctionSpaceType::dimRange;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      template< class GridFunction, std::enable_if_t< std::is_same< This, typename GridFunction::LocalFunctionType >::value, int > = 0 >
      SimpleLocalFunction ( const GridFunction &gridFunction )
        : localEvaluator_( gridFunction.localEvaluator() ), order_( gridFunction.order() )
      {}

      SimpleLocalFunction ( const EntityType &entity, LocalEvaluator localEvaluator, int order )
        : entity_( &entity ), localEvaluator_( std::move( localEvaluator ) ), order_( order )
      {}

      void init ( const EntityType &entity ) { entity_ = &entity; }

      template< class Point >
      void evaluate ( const Point &x, RangeType &value ) const
      {
        using Fem::coordinate;
        value = localEvaluator_( entity(), coordinate( x ) );
      }

      template< class Quadrature, class Values >
      void evaluateQuadrature ( const Quadrature &quadrature, Values &values ) const
      {
        for( const auto qp : quadrature )
          evaluate( qp, values[ qp.index() ] );
      }

      template< class Point >
      void jacobian ( const Point &x, JacobianRangeType &jacobian ) const
      {
        DUNE_THROW( NotImplemented, "SimpleLocalFunction::jacobian not implemented" );
      }

      template< class Quadrature, class Jacobians >
      void jacobianQuadrature ( const Quadrature &quadrature, Jacobians &jacobians ) const
      {
        for( const auto qp : quadrature )
          jacobian( qp, jacobians[ qp.index() ] );
      }

      template< class Point >
      void hessian ( const Point &x, HessianRangeType &hessian ) const
      {
        DUNE_THROW( NotImplemented, "SimpleLocalFunction::hessian not implemented" );
      }

      int order () const { return order_; }

      const EntityType &entity () const { assert( entity_ ); return *entity_; }

    private:
      const EntityType *entity_ = nullptr;
      LocalEvaluator localEvaluator_;
      int order_;
    };



    // IsLocalEvaluator
    // ----------------

    template< class LE, class E, class X >
    std::true_type __isLocalEvaluator ( const LE &, const E &, const X &, decltype( std::declval< LE >()( std::declval< E >(), std::declval< X >() ) ) * = nullptr );

    std::false_type __isLocalEvaluator ( ... );

    template< class GP, class LE >
    struct IsLocalEvaluator
      : public decltype( __isLocalEvaluator( std::declval< LE >(), std::declval< typename GP::template Codim< 0 >::EntityType >(), std::declval< typename GP::template Codim< 0 >::GeometryType::LocalCoordinate >() ) )
    {};



    // SimpleGridFunction
    // ------------------

    template< class GridPart, class LocalEvaluator >
    class SimpleGridFunction
      : public Fem::Function< typename SimpleLocalFunction< GridPart, LocalEvaluator >::FunctionSpaceType, SimpleGridFunction< GridPart, LocalEvaluator > >,
        public Fem::HasLocalFunction
    {
      typedef SimpleGridFunction< GridPart, LocalEvaluator > This;
      typedef Fem::Function< typename SimpleLocalFunction< GridPart, LocalEvaluator >::FunctionSpaceType, SimpleGridFunction< GridPart, LocalEvaluator > > Base;

      struct Space
      {
        Space(const GridPart &gridPart, int o)
          : gp_(gridPart), o_(o) {}
        int order() const
        {
          return o_;
        }
        const GridPart& gridPart() const
        {
          return gp_;
        }
        const GridPart &gp_;
        int o_;
      };
    public:
      typedef GridPart GridPartType;

      typedef SimpleLocalFunction< GridPart, LocalEvaluator > LocalFunctionType;

      typedef typename LocalFunctionType::EntityType EntityType;

      typedef typename Base::DomainType DomainType;
      typedef typename Base::RangeType RangeType;
      typedef typename Base::JacobianRangeType JacobianRangeType;

      SimpleGridFunction ( std::string name, const GridPartType &gridPart, LocalEvaluator localEvaluator, int order) // = std::numeric_limits< int >::max() )
        : name_( std::move( name ) ),
          gridPart_( gridPart ),
          localEvaluator_( std::move( localEvaluator ) ),
          order_( order )
      {}

      SimpleGridFunction ( const GridPartType &gridPart, LocalEvaluator localEvaluator, int order) // = std::numeric_limits< int >::max() )
        : name_( "unnamed-" + className< LocalEvaluator >() ),
          gridPart_( gridPart ),
          localEvaluator_( std::move( localEvaluator ) ),
          order_( order )
      {}

      LocalFunctionType localFunction ( const EntityType &entity ) const { return LocalFunctionType( entity, localEvaluator_, order_ ); }

      const std::string &name () const { return name_; }

      const GridPartType &gridPart () const { return gridPart_; }

      void evaluate ( const DomainType &x, RangeType &value ) const
      {
        Fem::EntitySearch< GridPartType, 0 > entitySearch( gridPart() );
        const EntityType entity = entitySearch( x );
        const auto geometry = entity.geometry();
        localFunction( entity ).evaluate( geometry.local( x ), value );
        value = RangeType(0);
      }

      void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const
      {
        Fem::EntitySearch< GridPartType, 0 > entitySearch( gridPart() );
        const EntityType entity = entitySearch( x );
        const auto geometry = entity.geometry();
        localFunction( entity ).jacobian( geometry.local( x ), jacobian );
        jacobian = JacobianRangeType(0);
      }

      const Space space() const
      {
        return Space(gridPart(), order_);
      }

      int order () const { return order_; }

      const LocalEvaluator &localEvaluator () const { return localEvaluator_; }

    protected:
      std::string name_;
      const GridPartType &gridPart_;
      LocalEvaluator localEvaluator_;
      int order_;
    };



    // IsEvaluator
    // -----------

    template< class E, class X >
    std::true_type __isEvaluator ( const E &, const X &, decltype( std::declval< E >()( std::declval< X >() ) ) * = nullptr );

    std::false_type __isEvaluator ( ... );

    template< class GP, class E >
    struct IsEvaluator
      : public decltype( __isEvaluator( std::declval< E >(), std::declval< typename GP::template Codim< 0 >::GeometryType::GlobalCoordinate >() ) )
    {};



    // LocalEvaluatorAdapter
    // ---------------------

    template< class Entity, class Evaluator >
    struct LocalEvaluatorAdapter
    {
      typedef typename Entity::Geometry::GlobalCoordinate GlobalCoordinate;
      typedef typename Entity::Geometry::LocalCoordinate LocalCoordinate;

      typedef decltype( std::declval< Evaluator >()( std::declval< GlobalCoordinate >() ) ) Value;

      LocalEvaluatorAdapter ( Evaluator evaluator ) : evaluator_( std::move( evaluator ) ) {}

      Value operator () ( const GlobalCoordinate &x ) const { return evaluator_( x ); }
      Value operator () ( const Entity &entity, const LocalCoordinate &x ) const { return evaluator_( entity.geometry().global( x ) ); }

    private:
      Evaluator evaluator_;
    };



    // SimpleGlobalGridFunction
    // ------------------------

    template< class GridPart, class Evaluator >
    class SimpleGlobalGridFunction
      : public SimpleGridFunction< GridPart, LocalEvaluatorAdapter< typename GridPart::template Codim< 0 >::EntityType, Evaluator > >
    {
      typedef SimpleGlobalGridFunction< GridPart, Evaluator > This;
      typedef SimpleGridFunction< GridPart, LocalEvaluatorAdapter< typename GridPart::template Codim< 0 >::EntityType, Evaluator > > Base;

    public:
      typedef typename Base::GridPartType GridPartType;

      typedef typename Base::DomainType DomainType;
      typedef typename Base::RangeType RangeType;

      SimpleGlobalGridFunction ( std::string name, const GridPartType &gridPart, Evaluator evaluator, int order = std::numeric_limits< int >::max() )
        : Base( std::move( name ), gridPart, std::move( evaluator ), order )
      {}

      SimpleGlobalGridFunction ( const GridPartType &gridPart, Evaluator evaluator, int order = std::numeric_limits< int >::max() )
        : Base( "unnamed-" + className< Evaluator >(), gridPart, std::move( evaluator ), order )
      {}

      void evaluate ( const DomainType &x, RangeType &value ) const { value = localEvaluator_( x ); }

    protected:
      using Base::localEvaluator_;
    };



    // simpleGridFunction
    // ------------------

    template< class GridPart, class LocalEvaluator >
    inline static std::enable_if_t< IsLocalEvaluator< GridPart, LocalEvaluator >::value, SimpleGridFunction< GridPart, LocalEvaluator > >
    simpleGridFunction ( std::string name, const GridPart &gridPart, LocalEvaluator localEvaluator, int order)
    {
      return SimpleGridFunction< GridPart, LocalEvaluator >( std::move( name ), gridPart, std::move( localEvaluator ), order );
    }

    template< class GridPart, class LocalEvaluator >
    inline static std::enable_if_t< IsLocalEvaluator< GridPart, LocalEvaluator >::value, SimpleGridFunction< GridPart, LocalEvaluator > >
    simpleGridFunction ( const GridPart &gridPart, LocalEvaluator localEvaluator, int order)
    {
      return SimpleGridFunction< GridPart, LocalEvaluator >( gridPart, std::move( localEvaluator ), order );
    }


    // simpleGridFunction
    // ------------------

    template< class GridPart, class Evaluator >
    inline static std::enable_if_t< IsEvaluator< GridPart, Evaluator >::value, SimpleGlobalGridFunction< GridPart, Evaluator > >
    simpleGridFunction ( std::string name, const GridPart &gridPart, Evaluator evaluator, int order)
    {
      return SimpleGlobalGridFunction< GridPart, Evaluator >( std::move( name ), gridPart, std::move( evaluator ), order );
    }

    template< class GridPart, class Evaluator >
    inline static std::enable_if_t< IsEvaluator< GridPart, Evaluator >::value, SimpleGlobalGridFunction< GridPart, Evaluator > >
    simpleGridFunction ( const GridPart &gridPart, Evaluator evaluator, int order)
    {
      return SimpleGlobalGridFunction< GridPart, Evaluator >( gridPart, std::move( evaluator ), order );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_FUNCTION_SIMPLEGRIDFUNCTION_HH
