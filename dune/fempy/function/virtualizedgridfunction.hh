#ifndef DUNE_FEMPY_FUNCTION_VIRTUALIZEDGRIDFUNCTION_HH
#define DUNE_FEMPY_FUNCTION_VIRTUALIZEDGRIDFUNCTION_HH

#include <functional>
#include <type_traits>
#include <utility>

#include <dune/common/visibility.hh>
#include <dune/common/ftraits.hh>

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/localfunction/bindable.hh>

#include <dune/fem/quadrature/elementquadrature.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fempy/quadrature/cachingpoint.hh>
#include <dune/fempy/quadrature/elementpoint.hh>
#include <dune/fempy/quadrature/fempyquadratures.hh>
#include <dune/fem/common/intersectionside.hh>

namespace Dune
{

  namespace FemPy
  {

    // VirtualizedGridFunction
    // ------------------------

    template< class GridPart, class Value >
    class VirtualizedGridFunction
      : public Fem::BindableGridFunctionWithSpace<GridPart,Dim<Value::dimension>>
    {
      typedef VirtualizedGridFunction< GridPart, Value > This;
      typedef Fem::BindableGridFunctionWithSpace<GridPart,Dim<Value::dimension>> Base;

    public:
      typedef typename GridPart::GridViewType GridView;
      typedef typename GridPart::template Codim< 0 >::EntityType EntityType;
      typedef typename GridPart::IntersectionType IntersectionType;

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
      template <int codim>
      using CachingPoint = FemPy::CachingPoint< GridPart, LocalCoordinateType, codim >;
      template <int codim>
      using ElementPoint = FemPy::ElementPoint< GridPart, LocalCoordinateType, codim >;

      // typically used quadratures for efficient evaluation of basis functions
      typedef Dune::Fem::CachingQuadrature< GridPart, 0 > ElementQuadratureType ;
      typedef Dune::Fem::CachingQuadrature< GridPart, 1 > FaceQuadratureType ;
      typedef Dune::Fem::CachingQuadrature< GridPart, 0, Dune::FemPy::FempyQuadratureTraits > ElementFemQuadratureType;
      typedef Dune::Fem::CachingQuadrature< GridPart, 1, Dune::FemPy::FempyQuadratureTraits > FaceFemQuadratureType;

      typedef std::vector< RangeType >          RangeValueVectorType;
      typedef std::vector< JacobianRangeType >  JacobianRangeValueVectorType;
      typedef std::vector< HessianRangeType >   HessianRangeValueVectorType;

      template< class QP >
      static Fem::QuadraturePointWrapper< QP > asQP ( const QP &qp )
      {
        return static_cast< Fem::QuadraturePointWrapper< QP > >( qp );
      }

      struct Interface
      {
        virtual ~Interface () = default;
        virtual Interface *clone () const = 0;

        virtual void evaluate ( const LocalCoordinateType &x, RangeType &value ) const = 0;
        virtual void evaluate ( const CachingPoint<0> &x, RangeType &value ) const = 0;
        virtual void evaluate ( const CachingPoint<1> &x, RangeType &value ) const = 0;
        virtual void evaluate ( const ElementPoint<0> &x, RangeType &value ) const = 0;
        virtual void evaluate ( const ElementPoint<1> &x, RangeType &value ) const = 0;
        virtual void jacobian ( const LocalCoordinateType &x, JacobianRangeType &jacobian ) const = 0;
        virtual void jacobian ( const CachingPoint<0> &x, JacobianRangeType &jacobian ) const = 0;
        virtual void jacobian ( const CachingPoint<1> &x, JacobianRangeType &jacobian ) const = 0;
        virtual void jacobian ( const ElementPoint<0> &x, JacobianRangeType &jacobian ) const = 0;
        virtual void jacobian ( const ElementPoint<1> &x, JacobianRangeType &jacobian ) const = 0;
        virtual void hessian ( const LocalCoordinateType &x, HessianRangeType &hessian ) const = 0;
        virtual void hessian ( const CachingPoint<0> &x, HessianRangeType &hessian ) const = 0;
        virtual void hessian ( const CachingPoint<1> &x, HessianRangeType &hessian ) const = 0;
        virtual void hessian ( const ElementPoint<0> &x, HessianRangeType &hessian ) const = 0;
        virtual void hessian ( const ElementPoint<1> &x, HessianRangeType &hessian ) const = 0;

        virtual void evaluateQuadrature( const ElementQuadratureType& quad, RangeValueVectorType& values ) const = 0 ;
        virtual void evaluateQuadrature( const FaceQuadratureType& quad, RangeValueVectorType& values ) const = 0 ;
        virtual void evaluateQuadrature( const ElementFemQuadratureType& quad, RangeValueVectorType& values ) const = 0 ;
        virtual void evaluateQuadrature( const FaceFemQuadratureType& quad, RangeValueVectorType& values ) const = 0 ;
        virtual void jacobianQuadrature( const ElementQuadratureType& quad, JacobianRangeValueVectorType& values ) const = 0 ;
        virtual void jacobianQuadrature( const FaceQuadratureType& quad, JacobianRangeValueVectorType& values ) const = 0 ;
        virtual void jacobianQuadrature( const ElementFemQuadratureType& quad, JacobianRangeValueVectorType& values ) const = 0 ;
        virtual void jacobianQuadrature( const FaceFemQuadratureType& quad, JacobianRangeValueVectorType& values ) const = 0 ;
        virtual void hessianQuadrature ( const ElementQuadratureType& quad, HessianRangeValueVectorType& values ) const = 0 ;
        virtual void hessianQuadrature ( const FaceQuadratureType& quad, HessianRangeValueVectorType& values ) const = 0 ;
        virtual void hessianQuadrature ( const ElementFemQuadratureType& quad, HessianRangeValueVectorType& values ) const = 0 ;
        virtual void hessianQuadrature ( const FaceFemQuadratureType& quad, HessianRangeValueVectorType& values ) const = 0 ;

        virtual void bind(const EntityType &entity) = 0;
        virtual void bind(const IntersectionType &intersection, Fem::IntersectionSide side) = 0;
        virtual void unbind() = 0;
      };

      template< class GF >
      struct DUNE_PRIVATE Implementation final
        : public Interface
      {
        typedef typename std::remove_reference<GF>::type GFType;

        Implementation ( const GF& gf ) :
          impl_( gf )
        {
        }

        virtual Interface *clone () const override { return new Implementation( *this ); }

        virtual void evaluate ( const LocalCoordinateType &x, RangeType &value ) const override { impl().evaluate( x, value ); }
        virtual void evaluate ( const CachingPoint<0> &x, RangeType &value ) const override { impl().evaluate( asQP( x ), value ); }
        virtual void evaluate ( const CachingPoint<1> &x, RangeType &value ) const override { impl().evaluate( asQP( x ), value ); }
        virtual void evaluate ( const ElementPoint<0> &x, RangeType &value ) const override { impl().evaluate( asQP( x ), value ); }
        virtual void evaluate ( const ElementPoint<1> &x, RangeType &value ) const override { impl().evaluate( asQP( x ), value ); }
        virtual void jacobian ( const LocalCoordinateType &x, JacobianRangeType &jacobian ) const override { impl().jacobian( x, jacobian ); }
        virtual void jacobian ( const CachingPoint<0> &x, JacobianRangeType &jacobian ) const override { impl().jacobian( asQP( x ), jacobian ); }
        virtual void jacobian ( const CachingPoint<1> &x, JacobianRangeType &jacobian ) const override { impl().jacobian( asQP( x ), jacobian ); }
        virtual void jacobian ( const ElementPoint<0> &x, JacobianRangeType &jacobian ) const override { impl().jacobian( asQP( x ), jacobian ); }
        virtual void jacobian ( const ElementPoint<1> &x, JacobianRangeType &jacobian ) const override { impl().jacobian( asQP( x ), jacobian ); }
        virtual void hessian ( const LocalCoordinateType &x, HessianRangeType &hessian ) const override { impl().hessian( x, hessian ); }
        virtual void hessian ( const CachingPoint<0> &x, HessianRangeType &hessian ) const override { impl().hessian( asQP( x ), hessian ); }
        virtual void hessian ( const CachingPoint<1> &x, HessianRangeType &hessian ) const override { impl().hessian( asQP( x ), hessian ); }
        virtual void hessian ( const ElementPoint<0> &x, HessianRangeType &hessian ) const override { impl().hessian( asQP( x ), hessian ); }
        virtual void hessian ( const ElementPoint<1> &x, HessianRangeType &hessian ) const override { impl().hessian( asQP( x ), hessian ); }
        virtual void evaluateQuadrature( const ElementQuadratureType& quad, RangeValueVectorType& values ) const override { impl().evaluateQuadrature( quad, values ); }
        virtual void evaluateQuadrature( const FaceQuadratureType& quad, RangeValueVectorType& values ) const override { impl().evaluateQuadrature( quad, values ); }
        virtual void evaluateQuadrature( const ElementFemQuadratureType& quad, RangeValueVectorType& values ) const override { impl().evaluateQuadrature( quad, values ); }
        virtual void evaluateQuadrature( const FaceFemQuadratureType& quad, RangeValueVectorType& values ) const override { impl().evaluateQuadrature( quad, values ); }
        virtual void jacobianQuadrature( const ElementQuadratureType& quad, JacobianRangeValueVectorType& values ) const override { impl().jacobianQuadrature( quad, values ); }
        virtual void jacobianQuadrature( const FaceQuadratureType& quad, JacobianRangeValueVectorType& values ) const override { impl().jacobianQuadrature( quad, values ); }
        virtual void jacobianQuadrature( const ElementFemQuadratureType& quad, JacobianRangeValueVectorType& values ) const override { impl().jacobianQuadrature( quad, values ); }
        virtual void jacobianQuadrature( const FaceFemQuadratureType& quad, JacobianRangeValueVectorType& values ) const override { impl().jacobianQuadrature( quad, values ); }
        virtual void hessianQuadrature ( const ElementQuadratureType& quad, HessianRangeValueVectorType& values ) const override { impl().hessianQuadrature( quad, values ); }
        virtual void hessianQuadrature ( const FaceQuadratureType& quad, HessianRangeValueVectorType& values ) const override { impl().hessianQuadrature( quad, values ); }
        virtual void hessianQuadrature ( const ElementFemQuadratureType& quad, HessianRangeValueVectorType& values ) const override { impl().hessianQuadrature( quad, values ); }
        virtual void hessianQuadrature ( const FaceFemQuadratureType& quad, HessianRangeValueVectorType& values ) const override { impl().hessianQuadrature( quad, values ); }

        virtual void bind(const EntityType &entity) override { impl_.bind(entity); }
        void bind(const IntersectionType &intersection, Fem::IntersectionSide side) override
        { impl_.bind(intersection, side); }
        virtual void unbind() override { impl_.unbind(); }
      private:
        const auto &impl () const { return impl_; }
        auto &impl () { return impl_; }

        Fem::ConstLocalFunction< GFType > impl_;
      };

    public:
      template <class Impl>
      VirtualizedGridFunction (const Impl &gf)
      : Base(gf.gridPart(),gf.name(),gf.order())
      , impl_( new Implementation< Impl >( gf ) )
      {}

      VirtualizedGridFunction ( const VirtualizedGridFunction &other )
      : Base(other.gridPart(),other.name(),other.order())
      , impl_( other ? other.impl_->clone() : nullptr )
      {}

      VirtualizedGridFunction ( VirtualizedGridFunction && ) = default;

      VirtualizedGridFunction &operator= ( const VirtualizedGridFunction &other ) { impl_.reset( other ? other.impl_->clone() : nullptr ); }
      VirtualizedGridFunction &operator= ( VirtualizedGridFunction && ) = default;

      explicit operator bool () const { return static_cast< bool >( impl_ ); }

      void bind(const EntityType &entity)
      {
        Base::bind(entity);
        impl_->bind(entity);
      }
      void bind(const IntersectionType &intersection, Fem::IntersectionSide side)
      {
        Base::bind(intersection, side);
        impl_->bind(intersection, side);
      }
      void unbind()
      {
        Base::unbind();
        impl_->unbind();
      }

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
        impl_->evaluate( CachingPoint<Quadrature::codimension>( x ), value );
      }
      template< class Quadrature >
      std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      evaluate ( const Fem::QuadraturePointWrapper< Quadrature > &x, RangeType &value ) const
      {
        impl_->evaluate( ElementPoint<Quadrature::codimension>( x ), value );
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
        impl_->jacobian( CachingPoint<Quadrature::codimension>( x ), jacobian );
      }
      template< class Quadrature >
      std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      jacobian ( const Fem::QuadraturePointWrapper< Quadrature > &x, JacobianRangeType &jacobian ) const
      {
        impl_->jacobian( ElementPoint<Quadrature::codimension>( x ), jacobian );
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
        impl_->hessian( CachingPoint<Quadrature::codimension>( x ), hessian );
      }
      template< class Quadrature >
      std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      hessian ( const Fem::QuadraturePointWrapper< Quadrature > &x, HessianRangeType &hessian ) const
      {
        impl_->hessian( ElementPoint<Quadrature::codimension>( x ), hessian );
      }

      template< class Quadrature, class ... Vectors >
      void evaluateQuadrature ( const Quadrature &quad, Vectors & ... values ) const
      {
        static_assert( sizeof...( Vectors ) > 0, "evaluateQuadrature needs to be called with at least one vector." );
        std::ignore = std::make_tuple( ( evaluateSingleQuadrature( quad, values ), 1 ) ... );
      }

      void evaluateSingleQuadrature ( const ElementQuadratureType &quadrature, RangeValueVectorType &values ) const
      { impl_->evaluateQuadrature( quadrature, values ); }
      void evaluateSingleQuadrature ( const FaceQuadratureType &quadrature, RangeValueVectorType &values ) const
      { impl_->evaluateQuadrature( quadrature, values ); }
      void evaluateSingleQuadrature ( const ElementFemQuadratureType &quadrature, RangeValueVectorType &values ) const
      { impl_->evaluateQuadrature( quadrature, values ); }
      void evaluateSingleQuadrature ( const FaceFemQuadratureType &quadrature, RangeValueVectorType &values ) const
      { impl_->evaluateQuadrature( quadrature, values ); }
      template< class Quadrature, class Vector >
      auto evaluateSingleQuadrature ( const Quadrature &quad, Vector &v ) const
      -> std::enable_if_t< std::is_same< std::decay_t< decltype(v[ 0 ]) >, RangeType >::value >
      {
        for( const auto qp : quad )
          evaluate( qp, v[ qp.index() ] );
      }

      void jacobianQuadrature ( const ElementQuadratureType &quadrature, JacobianRangeValueVectorType &values ) const
      { impl_->jacobianQuadrature( quadrature, values ); }
      void jacobianQuadrature ( const FaceQuadratureType &quadrature, JacobianRangeValueVectorType &values ) const
      { impl_->jacobianQuadrature( quadrature, values ); }
      void jacobianQuadrature ( const ElementFemQuadratureType &quadrature, JacobianRangeValueVectorType &values ) const
      { impl_->jacobianQuadrature( quadrature, values ); }
      void jacobianQuadrature ( const FaceFemQuadratureType &quadrature, JacobianRangeValueVectorType &values ) const
      { impl_->jacobianQuadrature( quadrature, values ); }

      template< class Quadrature, class Jacobians >
      void jacobianQuadrature ( const Quadrature &quadrature, Jacobians &jacobians ) const
      {
        for( const auto qp : quadrature )
          jacobian( qp, jacobians[ qp.index() ] );
      }
      template< class Quadrature, class Vector >
      auto evaluateSingleQuadrature ( const Quadrature &quad, Vector &v ) const
      -> std::enable_if_t< std::is_same< std::decay_t< decltype(v[ 0 ]) >, JacobianRangeType >::value >
      { jacobianQuadrature(quad,v); }

      void hessianQuadrature ( const ElementQuadratureType &quadrature, HessianRangeValueVectorType &values ) const
      { impl_->hessianQuadrature( quadrature, values ); }
      void hessianQuadrature ( const FaceQuadratureType &quadrature, HessianRangeValueVectorType &values ) const
      { impl_->hessianQuadrature( quadrature, values ); }
      void hessianQuadrature ( const ElementFemQuadratureType &quadrature, HessianRangeValueVectorType &values ) const
      { impl_->hessianQuadrature( quadrature, values ); }
      void hessianQuadrature ( const FaceFemQuadratureType &quadrature, HessianRangeValueVectorType &values ) const
      { impl_->hessianQuadrature( quadrature, values ); }

      template< class Quadrature, class Hessians >
      void hessianQuadrature ( const Quadrature &quadrature, Hessians &hessians ) const
      {
        for( const auto qp : quadrature )
          hessian( qp, hessians[ qp.index() ] );
      }
      template< class Quadrature, class Vector >
      auto evaluateSingleQuadrature ( const Quadrature &quad, Vector &v ) const
      -> std::enable_if_t< std::is_same< std::decay_t< decltype(v[ 0 ]) >, HessianRangeType >::value >
      { hessianQuadrature(quad,v); }

    private:
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
