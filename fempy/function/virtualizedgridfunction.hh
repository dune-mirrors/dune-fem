#ifndef DUNE_FEMPY_FUNCTION_VIRTUALIZEDGRIDFUNCTION_HH
#define DUNE_FEMPY_FUNCTION_VIRTUALIZEDGRIDFUNCTION_HH

#include <functional>
#include <type_traits>
#include <utility>

#include <dune/common/ftraits.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/space/common/functionspace.hh>

#include <dune/fempy/quadrature/cachingpoint.hh>

namespace Dune
{

  namespace FemPy
  {

    // VirtualizedLocalFunction
    // ------------------------

    template< class GridPart, class Value >
    class VirtualizedLocalFunction
    {
      typedef VirtualizedLocalFunction< GridPart, Value > This;

    public:
      typedef typename GridPart::template Codim< 0 >::EntityType EntityType;

      typedef typename EntityType::Geometry::LocalCoordinate LocalCoordinateType;
      typedef typename EntityType::Geometry::GlobalCoordinate GlobalCoordinateType;

      typedef Fem::FunctionSpace< typename GridPart::ctype, typename FieldTraits< Value >::field_type, GridPart::dimensionworld, Value::dimension > FunctionSpaceType;

      static const int dimDomain = FunctionSpaceType::dimDomain;
      static const int dimRange = FunctionSpaceType::dimRange;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    private:
      typedef Fem::QuadraturePointWrapper< CachingPoint< LocalCoordinateType > > QuadraturePoint;

      template< class GF >
      static constexpr bool isGridFunction ()
      {
        using std::cref;
        return std::is_base_of< Fem::HasLocalFunction, std::decay_t< decltype( cref( std::declval< const GF & >() ).get() ) > >::value;
      }

      template< class GF >
      static std::true_type canCreateLocalFunctionHelp ( const GF &, std::decay_t< decltype( std::declval< GF >().localFunction() ) > * = nullptr );
      static std::false_type canCreateLocalFunctionHelp ( ... );

      template< class GF >
      static constexpr bool canCreateLocalFunction ()
      {
        return decltype( canCreateLocalFunctionHelp ( std::declval< const GF & >() ) )::value;
      }

      struct Interface
      {
        virtual ~Interface () = default;
        virtual Interface *clone () const = 0;

        virtual void init ( const EntityType &entity ) = 0;
        virtual void evaluate ( const LocalCoordinateType &x, RangeType &value ) const = 0;
        virtual void evaluate ( const QuadraturePoint &x, RangeType &value ) const = 0;
        virtual void jacobian ( const LocalCoordinateType &x, JacobianRangeType &jacobian ) const = 0;
        virtual void jacobian ( const QuadraturePoint &x, JacobianRangeType &jacobian ) const = 0;
        virtual int order () const = 0;
        virtual const EntityType &entity () const = 0;
      };

      template< class Impl >
      struct Implementation final
        : public Interface
      {
        Implementation ( Impl impl ) : impl_( std::move( impl ) ) {}
        virtual Interface *clone () const override { return new Implementation( *this ); }

        virtual void init ( const EntityType &entity ) override { impl().init( entity ); }

        virtual void evaluate ( const LocalCoordinateType &x, RangeType &value ) const override { impl().evaluate( x, value ); }
        virtual void evaluate ( const QuadraturePoint &x, RangeType &value ) const override { impl().evaluate( x, value ); }
        virtual void jacobian ( const LocalCoordinateType &x, JacobianRangeType &jacobian ) const override { impl().jacobian( x, jacobian ); }
        virtual void jacobian ( const QuadraturePoint &x, JacobianRangeType &jacobian ) const override { impl().jacobian( x, jacobian ); }
        virtual int order () const override { return impl().order(); }
        virtual const EntityType &entity () const override { return impl().entity(); }

      private:
        const auto &impl () const { using std::cref; return cref( impl_ ).get(); }
        auto &impl () { using std::ref; return ref( impl_ ).get(); }

        Impl impl_;
      };

    public:
      template< class Impl, std::enable_if_t< !isGridFunction< Impl >() && !std::is_base_of< VirtualizedLocalFunction, Impl >::value, int > = 0 >
      VirtualizedLocalFunction ( Impl impl )
        : impl_( new Implementation< Impl >( std::move( impl ) ) )
      {}

      template< class GF, std::enable_if_t< isGridFunction< GF >() && canCreateLocalFunction< GF >(), int > = 0 >
      VirtualizedLocalFunction ( const GF &gf )
        : VirtualizedLocalFunction( gf.localFunction() )
      {}

      template< class GF, std::enable_if_t< isGridFunction< GF >() && !canCreateLocalFunction< GF >(), int > = 0 >
      VirtualizedLocalFunction ( const GF &gf )
        : VirtualizedLocalFunction( typename GF::LocalFunctionType( gf ) )
      {}

      VirtualizedLocalFunction ( const VirtualizedLocalFunction &other ) : impl_( other.impl_->clone() ) {}
      VirtualizedLocalFunction ( VirtualizedLocalFunction && ) = default;

      void init ( const EntityType &entity ) { impl_->init( entity ); }

      template< class Point >
      void evaluate ( const Point &x, RangeType &value ) const
      {
        using Fem::coordinate;
        impl_->evaluate( coordinate( x ), value );
      }

      template< class Quadrature >
      std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      evaluate ( const Fem::QuadraturePointWrapper< Quadrature > &x, RangeType &value ) const
      {
        // note: QuadraturePoint will store a reference to QuadratureCacheIdentifier.
        //       This object must be kept alive until the evaluation returns!
        CachingPoint< LocalCoordinateType > y( x );
        impl_->evaluate( static_cast< QuadraturePoint >( y ), value );
      }

      template< class Quadrature >
      std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      evaluate ( const Fem::QuadraturePointWrapper< Quadrature > &x, RangeType &value ) const
      {
        impl_->evaluate( x.position(), value );
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
        using Fem::coordinate;
        impl_->jacobian( coordinate( x ), jacobian );
      }

      template< class Quadrature >
      std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      jacobian ( const Fem::QuadraturePointWrapper< Quadrature > &x, JacobianRangeType &jacobian ) const
      {
        // note: QuadraturePoint will store a reference to QuadratureCacheIdentifier.
        //       This object must be kept alive until the evaluation returns!
        CachingPoint< LocalCoordinateType > y( x );
        impl_->jacobian( static_cast< QuadraturePoint >( y ), jacobian );
      }

      template< class Quadrature >
      std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      jacobian ( const Fem::QuadraturePointWrapper< Quadrature > &x, JacobianRangeType &jacobian ) const
      {
        impl_->jacobian( x.position(), jacobian );
      }

      template< class Quadrature, class Jacobians >
      void jacobianQuadrature ( const Quadrature &quadrature, Jacobians &jacobians ) const
      {
        for( const auto qp : quadrature )
          jacobian( qp, jacobians[ qp.index() ] );
      }

      int order () const { return impl_->order(); }

      const EntityType &entity () const { return impl_->entity(); }

    private:
      std::unique_ptr< Interface > impl_;
    };



    // VirtualizedGridFunction
    // -----------------------

    template< class GridPart, class Value >
    class VirtualizedGridFunction
      : public Fem::Function< typename VirtualizedLocalFunction< GridPart, Value >::FunctionSpaceType, VirtualizedGridFunction< GridPart, Value > >,
        public Fem::HasLocalFunction
    {
      typedef VirtualizedGridFunction< GridPart, Value > This;
      typedef Fem::Function< typename VirtualizedLocalFunction< GridPart, Value >::FunctionSpaceType, VirtualizedGridFunction< GridPart, Value > > Base;

    public:
      typedef GridPart GridPartType;

      typedef VirtualizedLocalFunction< GridPart, Value > LocalFunctionType;

      typedef typename LocalFunctionType::EntityType EntityType;

      typedef typename Base::DomainType DomainType;
      typedef typename Base::RangeType RangeType;
      typedef typename Base::JacobianRangeType JacobianRangeType;

    private:
      template< class GF >
      static constexpr bool isGridFunction ()
      {
        using std::cref;
        return std::is_base_of< Fem::HasLocalFunction, std::decay_t< decltype( cref( std::declval< const GF & >() ).get() ) > >::value;
      }

      struct Interface
      {
        virtual ~Interface () = default;
        virtual Interface *clone () const = 0;

        virtual LocalFunctionType localFunction () const = 0;
        virtual LocalFunctionType localFunction ( const EntityType &entity ) const = 0;
        virtual std::string name () const = 0;
        virtual const GridPartType &gridPart () const = 0;
        virtual void evaluate ( const DomainType &x, RangeType &value ) const = 0;
        virtual void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const = 0;
      };

      template< class Impl >
      struct Implementation final
        : public Interface
      {
        Implementation ( Impl impl ) : impl_( std::move( impl ) ) {}
        virtual Interface *clone () const override { return new Implementation( *this ); }

        virtual LocalFunctionType localFunction () const override { return std::decay_t< decltype( impl().localFunction( std::declval< const EntityType & >() ) ) >( impl() ); }
        virtual LocalFunctionType localFunction ( const EntityType &entity ) const override { return impl().localFunction( entity ); }
        virtual std::string name () const override { return impl().name(); }
        virtual const GridPartType &gridPart () const override { return impl().gridPart(); }
        virtual void evaluate ( const DomainType &x, RangeType &value ) const override { impl().evaluate( x, value ); }
        virtual void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const override { impl().jacobian( x, jacobian ); }

      private:
        auto &impl () const { using std::cref; return cref( impl_ ).get(); }

        Impl impl_;
      };

    public:
      template< class Impl, std::enable_if_t< isGridFunction< Impl >() && !std::is_base_of< VirtualizedGridFunction, Impl >::value, int > = 0 >
      VirtualizedGridFunction ( Impl impl )
        : impl_( new Implementation< Impl >( std::move( impl ) ) )
      {}

      VirtualizedGridFunction ( const VirtualizedGridFunction &other ) : impl_( other.impl_->clone() ) {}
      VirtualizedGridFunction ( VirtualizedGridFunction && ) = default;

      LocalFunctionType localFunction () const { return impl_->localFunction(); }
      LocalFunctionType localFunction ( const EntityType &entity ) const { return impl_->localFunction( entity ); }

      std::string name () const { return impl_->name(); }

      const GridPartType &gridPart () const { return impl_->gridPart(); }

      void evaluate ( const DomainType &x, RangeType &value ) const { return impl_->evaluate( x, value ); }
      void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const { return impl_->jacobian( x, jacobian ); }

    protected:
      std::unique_ptr< Interface > impl_;
    };



    // virtualizedGridFunction
    // -----------------------

    template< class GridFunction >
    inline static auto virtualizeGridFunction ( GridFunction gridFunction )
    {
      typedef std::decay_t< decltype( std::ref( std::declval< GridFunction >() ).get() ) > Impl;
      return VirtualizedGridFunction< typename Impl::GridPartType, typename Impl::RangeType >( std::move( gridFunction ) );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_FUNCTION_VIRTUALIZEDGRIDFUNCTION_HH
