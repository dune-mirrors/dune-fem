#ifndef DUNE_FEMPY_FUNCTION_VIRTUALIZEDGRIDFUNCTION_HH
#define DUNE_FEMPY_FUNCTION_VIRTUALIZEDGRIDFUNCTION_HH

#include <functional>
#include <type_traits>
#include <utility>

#include <dune/common/ftraits.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/space/common/functionspace.hh>

#include <dune/fempy/quadrature/cachingpoint.hh>
#include <dune/fempy/quadrature/elementpoint.hh>

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

      typedef typename FieldTraits< Value >::field_type RangeFieldType;
      typedef Fem::FunctionSpace< typename GridPart::ctype, RangeFieldType, GridPart::dimensionworld, Value::dimension > FunctionSpaceType;

      static const int dimDomain = FunctionSpaceType::dimDomain;
      static const int dimRange = FunctionSpaceType::dimRange;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

    private:
      typedef FemPy::CachingPoint< LocalCoordinateType > CachingPoint;
      typedef FemPy::ElementPoint< LocalCoordinateType > ElementPoint;

      // typically used quadratures for efficient evaluation of basis functions
      typedef Dune::Fem::CachingQuadrature< GridPart, 0 > ElementQuadratureType ;
      typedef Dune::Fem::CachingQuadrature< GridPart, 1 > FaceQuadratureType ;

      typedef std::vector< RangeType >          RangeValueVectorType;
      typedef std::vector< JacobianRangeType >  JacobianRangeValueVectorType;
      typedef std::vector< HessianRangeType >   HessianRangeValueVectorType;

      template< class QP >
      static Fem::QuadraturePointWrapper< QP > asQP ( const QP &qp )
      {
        return static_cast< Fem::QuadraturePointWrapper< QP > >( qp );
      }

      template< class GF >
      static constexpr bool isGridFunction ()
      {
        using std::cref;
        return std::is_base_of< Fem::HasLocalFunction, std::decay_t< decltype( cref( std::declval< const GF & >() ).get() ) > >::value;
      }

      template< class LF >
      static std::true_type isLocalFunctionHelp ( const LF &, decltype( std::declval< const LF & >().evaluate( std::declval< const typename LF::EntityType::Geometry::LocalCoordinate & >(), std::declval< typename LF::RangeType & >() ) ) * = nullptr );
      static std::false_type isLocalFunctionHelp ( ... );

      template< class LF >
      static constexpr bool isLocalFunction ()
      {
        return (!isGridFunction< LF >() && decltype( isLocalFunctionHelp( std::declval< const LF & >() ) )::value);
      }

      template< class GF >
      static std::true_type canCreateLocalFunctionHelp ( const GF &, std::decay_t< decltype( std::declval< GF >().localFunction() ) > * = nullptr );
      static std::false_type canCreateLocalFunctionHelp ( ... );

      template< class GF >
      static constexpr bool canCreateLocalFunction ()
      {
        return decltype( canCreateLocalFunctionHelp( std::declval< const GF & >() ) )::value;
      }

      struct Interface
      {
        virtual ~Interface () = default;
        virtual Interface *clone () const = 0;

        virtual void init ( const EntityType &entity ) = 0;
        virtual void evaluate ( const LocalCoordinateType &x, RangeType &value ) const = 0;
        virtual void evaluate ( const CachingPoint &x, RangeType &value ) const = 0;
        virtual void evaluate ( const ElementPoint &x, RangeType &value ) const = 0;
        virtual void jacobian ( const LocalCoordinateType &x, JacobianRangeType &jacobian ) const = 0;
        virtual void jacobian ( const CachingPoint &x, JacobianRangeType &jacobian ) const = 0;
        virtual void jacobian ( const ElementPoint &x, JacobianRangeType &jacobian ) const = 0;
        virtual void hessian ( const LocalCoordinateType &x, HessianRangeType &hessian ) const = 0;
        virtual void hessian ( const CachingPoint &x, HessianRangeType &hessian ) const = 0;
        virtual void hessian ( const ElementPoint &x, HessianRangeType &hessian ) const = 0;

        virtual void evaluateQuadrature( const ElementQuadratureType& quad, RangeValueVectorType& values ) const = 0 ;
        virtual void evaluateQuadrature( const FaceQuadratureType& quad, RangeValueVectorType& values ) const = 0 ;
        virtual void jacobianQuadrature( const ElementQuadratureType& quad, JacobianRangeValueVectorType& values ) const = 0 ;
        virtual void jacobianQuadrature( const FaceQuadratureType& quad, JacobianRangeValueVectorType& values ) const = 0 ;
        virtual void hessianQuadrature ( const ElementQuadratureType& quad, HessianRangeValueVectorType& values ) const = 0 ;
        virtual void hessianQuadrature ( const FaceQuadratureType& quad, HessianRangeValueVectorType& values ) const = 0 ;

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
        virtual void evaluate ( const CachingPoint &x, RangeType &value ) const override { impl().evaluate( asQP( x ), value ); }
        virtual void evaluate ( const ElementPoint &x, RangeType &value ) const override { impl().evaluate( asQP( x ), value ); }
        virtual void jacobian ( const LocalCoordinateType &x, JacobianRangeType &jacobian ) const override { impl().jacobian( x, jacobian ); }
        virtual void jacobian ( const CachingPoint &x, JacobianRangeType &jacobian ) const override { impl().jacobian( asQP( x ), jacobian ); }
        virtual void jacobian ( const ElementPoint &x, JacobianRangeType &jacobian ) const override { impl().jacobian( asQP( x ), jacobian ); }
        virtual void hessian ( const LocalCoordinateType &x, HessianRangeType &hessian ) const override { impl().hessian( x, hessian ); }
        virtual void hessian ( const CachingPoint &x, HessianRangeType &hessian ) const override { impl().hessian( asQP( x ), hessian ); }
        virtual void hessian ( const ElementPoint &x, HessianRangeType &hessian ) const override { impl().hessian( asQP( x ), hessian ); }
        virtual void evaluateQuadrature( const ElementQuadratureType& quad, RangeValueVectorType& values ) const { impl().evaluateQuadrature( quad, values ); }
        virtual void evaluateQuadrature( const FaceQuadratureType& quad, RangeValueVectorType& values ) const { impl().evaluateQuadrature( quad, values ); }
        virtual void jacobianQuadrature( const ElementQuadratureType& quad, JacobianRangeValueVectorType& values ) const { impl().jacobianQuadrature( quad, values ); }
        virtual void jacobianQuadrature( const FaceQuadratureType& quad, JacobianRangeValueVectorType& values ) const { impl().jacobianQuadrature( quad, values ); }
        virtual void hessianQuadrature ( const ElementQuadratureType& quad, HessianRangeValueVectorType& values ) const { impl().hessianQuadrature( quad, values ); }
        virtual void hessianQuadrature ( const FaceQuadratureType& quad, HessianRangeValueVectorType& values ) const { impl().hessianQuadrature( quad, values ); }

        virtual int order () const override { return impl().order(); }
        virtual const EntityType &entity () const override { return impl().entity(); }

      private:
        const auto &impl () const { using std::cref; return cref( impl_ ).get(); }
        auto &impl () { using std::ref; return ref( impl_ ).get(); }

        Impl impl_;
      };

    public:
      VirtualizedLocalFunction () = default;

      template< class LF, std::enable_if_t< isLocalFunction< LF >() && !std::is_base_of< VirtualizedLocalFunction, LF >::value, int > = 0 >
      VirtualizedLocalFunction ( LF lf )
        : impl_( new Implementation< LF >( std::move( lf ) ) )
      {}

      template< class GF, std::enable_if_t< isGridFunction< GF >() && canCreateLocalFunction< GF >(), int > = 0 >
      explicit VirtualizedLocalFunction ( const GF &gf )
        : VirtualizedLocalFunction( gf.localFunction() )
      {}

      template< class GF, std::enable_if_t< isGridFunction< GF >() && !canCreateLocalFunction< GF >(), int > = 0 >
      explicit VirtualizedLocalFunction ( const GF &gf )
        : VirtualizedLocalFunction( typename GF::LocalFunctionType( gf ) )
      {}

      VirtualizedLocalFunction ( const VirtualizedLocalFunction &other ) : impl_( other ? other.impl_->clone() : nullptr ) {}
      VirtualizedLocalFunction ( VirtualizedLocalFunction && ) = default;

      VirtualizedLocalFunction &operator= ( const VirtualizedLocalFunction &other ) { impl_.reset( other ? other.impl_->clone() : nullptr ); }
      VirtualizedLocalFunction &operator= ( VirtualizedLocalFunction && ) = default;

      explicit operator bool () const { return static_cast< bool >( impl_ ); }

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
        impl_->evaluate( CachingPoint( x ), value );
      }

      template< class Quadrature >
      std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      evaluate ( const Fem::QuadraturePointWrapper< Quadrature > &x, RangeType &value ) const
      {
        impl_->evaluate( ElementPoint( x ), value );
      }

      template< class Quadrature, class Values >
      void evaluateQuadrature ( const Quadrature &quadrature, Values &values ) const
      {
        for( const auto qp : quadrature )
          evaluate( qp, values[ qp.index() ] );
      }

      void evaluateQuadrature ( const ElementQuadratureType &quadrature, RangeValueVectorType &values ) const
      {
        std::cout << "Virtualized::evalQuad " << std::endl;
        impl_->evaluateQuadrature( quadrature, values );
      }

      void evaluateQuadrature ( const FaceQuadratureType &quadrature, RangeValueVectorType &values ) const
      {
        std::cout << "Virtualized::evalQuad " << std::endl;
        impl_->evaluateQuadrature( quadrature, values );
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
        impl_->jacobian( CachingPoint( x ), jacobian );
      }

      template< class Quadrature >
      std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      jacobian ( const Fem::QuadraturePointWrapper< Quadrature > &x, JacobianRangeType &jacobian ) const
      {
        impl_->jacobian( ElementPoint( x ), jacobian );
      }

      template< class Quadrature, class Jacobians >
      void jacobianQuadrature ( const Quadrature &quadrature, Jacobians &jacobians ) const
      {
        for( const auto qp : quadrature )
          jacobian( qp, jacobians[ qp.index() ] );
      }

      template< class Quadrature, class Hessians >
      void hessianQuadrature ( const Quadrature &quadrature, Hessians &hessians ) const
      {
        for( const auto qp : quadrature )
          hessian( qp, hessians[ qp.index() ] );
      }

      template< class Point >
      void hessian ( const Point &x, HessianRangeType &hessian ) const
      {
        using Fem::coordinate;
        impl_->hessian( coordinate( x ), hessian );
      }

      template< class Quadrature >
      std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      hessian ( const Fem::QuadraturePointWrapper< Quadrature > &x, HessianRangeType &hessian ) const
      {
        impl_->hessian( CachingPoint( x ), hessian );
      }

      template< class Quadrature >
      std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      hessian ( const Fem::QuadraturePointWrapper< Quadrature > &x, HessianRangeType &hessian ) const
      {
        impl_->hessian( ElementPoint( x ), hessian );
      }

      int order () const { return impl_->order(); }

      const EntityType &entity () const { assert( impl_ ); return impl_->entity(); }

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

      static const int dimRange = RangeType::dimension;

    private:
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
      const Space space() const
      {
        return Space(impl_->gridPart(), 5);
      }

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
