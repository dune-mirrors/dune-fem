#ifndef DUNE_FEMPY_FUNCTION_VIRTUALIZEDCONSTLOCALFUNCTION_HH
#define DUNE_FEMPY_FUNCTION_VIRTUALIZEDCONSTLOCALFUNCTION_HH

#include <dune/fempy/function/virtualizedgridfunction.hh>

namespace Dune
{

  namespace FemPy
  {

    // VirtualizedConstLocalFunctionImpl
    // -----------------------------

    template< class GridPart, class Value, bool wrapper = true >
    class VirtualizedConstLocalFunctionImpl
    {
      typedef VirtualizedConstLocalFunctionImpl< GridPart, Value, wrapper > This;
      typedef VirtualizedGridFunction< GridPart, Value > VGF;

    public:
      typedef typename GridPart::GridViewType GridView;
      typedef typename GridPart::template Codim< 0 >::EntityType EntityType;
      typedef typename EntityType::Geometry                      Geometry;
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

    public:
      typedef typename VGF::Interface Interface;

      // stores a references to the GF (copy = ! wrapper)
      template <class GF>
      using Implementation = typename VGF::Implementation< GF, !wrapper /* copy */>;

    public:
      //! constructor taking a ConstLocalFunction< DiscreteFunction >
      template < class GF >
      VirtualizedConstLocalFunctionImpl (const GF &gf)
      : impl_( new Implementation< GF >( gf ) )
      {}

      /*
      template < class... Args >
      VirtualizedConstLocalFunctionImpl (const Fem::ConstLocalFunction< Args... > &gf)
      : impl_( new Implementation< Fem::ConstLocalFunction< Args... > >( gf ) )
      {}
      */

      VirtualizedConstLocalFunctionImpl ( const VirtualizedConstLocalFunctionImpl &other )
      : impl_( other ? other.impl_->clone() : nullptr )
      {}

      VirtualizedConstLocalFunctionImpl ( VirtualizedConstLocalFunctionImpl && ) = default;

      VirtualizedConstLocalFunctionImpl &operator= ( const VirtualizedConstLocalFunctionImpl &other ) { impl_.reset( other ? other.impl_->clone() : nullptr ); }
      VirtualizedConstLocalFunctionImpl &operator= ( VirtualizedConstLocalFunctionImpl && ) = default;

      explicit operator bool () const { return static_cast< bool >( impl_ ); }

      const EntityType &entity () const { return impl().entity(); }
      const Geometry &geometry () const { return impl().geometry(); }

      void bind(const EntityType &entity)
      {
        impl().bind(entity);
      }
      void bind(const IntersectionType &intersection, Fem::IntersectionSide side)
      {
        impl().bind(intersection, side);
      }
      void unbind()
      {
        impl().unbind();
      }

      template< class Point >
      void evaluate ( const Point &x, RangeType &value ) const
      {
        using Fem::coordinate;
        impl().evaluate( coordinate( x ), value );
      }
      template< class Quadrature >
      std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      evaluate ( const Fem::QuadraturePointWrapper< Quadrature > &x, RangeType &value ) const
      {
        impl().evaluate( CachingPoint<Quadrature::codimension>( x ), value );
      }
      template< class Quadrature >
      std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      evaluate ( const Fem::QuadraturePointWrapper< Quadrature > &x, RangeType &value ) const
      {
        impl().evaluate( ElementPoint<Quadrature::codimension>( x ), value );
      }

      template< class Point >
      void jacobian ( const Point &x, JacobianRangeType &jacobian ) const
      {
        using Fem::coordinate;
        impl().jacobian( coordinate( x ), jacobian );
      }
      template< class Quadrature >
      std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      jacobian ( const Fem::QuadraturePointWrapper< Quadrature > &x, JacobianRangeType &jacobian ) const
      {
        impl().jacobian( CachingPoint<Quadrature::codimension>( x ), jacobian );
      }
      template< class Quadrature >
      std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      jacobian ( const Fem::QuadraturePointWrapper< Quadrature > &x, JacobianRangeType &jacobian ) const
      {
        impl().jacobian( ElementPoint<Quadrature::codimension>( x ), jacobian );
      }

      template< class Point >
      void hessian ( const Point &x, HessianRangeType &hessian ) const
      {
        using Fem::coordinate;
        impl().hessian( coordinate( x ), hessian );
      }
      template< class Quadrature >
      std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      hessian ( const Fem::QuadraturePointWrapper< Quadrature > &x, HessianRangeType &hessian ) const
      {
        impl().hessian( CachingPoint<Quadrature::codimension>( x ), hessian );
      }
      template< class Quadrature >
      std::enable_if_t< !std::is_convertible< Quadrature, Fem::CachingInterface >::value >
      hessian ( const Fem::QuadraturePointWrapper< Quadrature > &x, HessianRangeType &hessian ) const
      {
        impl().hessian( ElementPoint<Quadrature::codimension>( x ), hessian );
      }

      template< class Quadrature, class ... Vectors >
      void evaluateQuadrature ( const Quadrature &quad, Vectors & ... values ) const
      {
        static_assert( sizeof...( Vectors ) > 0, "evaluateQuadrature needs to be called with at least one vector." );
        std::ignore = std::make_tuple( ( evaluateSingleQuadrature( quad, values ), 1 ) ... );
      }

      void evaluateSingleQuadrature ( const ElementQuadratureType &quadrature, RangeValueVectorType &values ) const
      { impl().evaluateQuadrature( quadrature, values ); }
      void evaluateSingleQuadrature ( const FaceQuadratureType &quadrature, RangeValueVectorType &values ) const
      { impl().evaluateQuadrature( quadrature, values ); }
      void evaluateSingleQuadrature ( const ElementFemQuadratureType &quadrature, RangeValueVectorType &values ) const
      { impl().evaluateQuadrature( quadrature, values ); }
      void evaluateSingleQuadrature ( const FaceFemQuadratureType &quadrature, RangeValueVectorType &values ) const
      { impl().evaluateQuadrature( quadrature, values ); }
      template< class Quadrature, class Vector >
      auto evaluateSingleQuadrature ( const Quadrature &quad, Vector &v ) const
      -> std::enable_if_t< std::is_same< std::decay_t< decltype(v[ 0 ]) >, RangeType >::value >
      {
        for( const auto qp : quad )
          evaluate( qp, v[ qp.index() ] );
      }

      void jacobianQuadrature ( const ElementQuadratureType &quadrature, JacobianRangeValueVectorType &values ) const
      { impl().jacobianQuadrature( quadrature, values ); }
      void jacobianQuadrature ( const FaceQuadratureType &quadrature, JacobianRangeValueVectorType &values ) const
      { impl().jacobianQuadrature( quadrature, values ); }
      void jacobianQuadrature ( const ElementFemQuadratureType &quadrature, JacobianRangeValueVectorType &values ) const
      { impl().jacobianQuadrature( quadrature, values ); }
      void jacobianQuadrature ( const FaceFemQuadratureType &quadrature, JacobianRangeValueVectorType &values ) const
      { impl().jacobianQuadrature( quadrature, values ); }

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
      { impl().hessianQuadrature( quadrature, values ); }
      void hessianQuadrature ( const FaceQuadratureType &quadrature, HessianRangeValueVectorType &values ) const
      { impl().hessianQuadrature( quadrature, values ); }
      void hessianQuadrature ( const ElementFemQuadratureType &quadrature, HessianRangeValueVectorType &values ) const
      { impl().hessianQuadrature( quadrature, values ); }
      void hessianQuadrature ( const FaceFemQuadratureType &quadrature, HessianRangeValueVectorType &values ) const
      { impl().hessianQuadrature( quadrature, values ); }

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

    protected:
      const Interface& impl() const { assert( impl_ ); return *impl_; }
      Interface& impl() { assert( impl_ ); return *impl_; }

      std::unique_ptr< Interface > impl_;
    };

    //! VirtualizedConstLocalFunctionWrapper stores the wrapped ConstLocalFunction as a reference
    template< class GridPart, class Value >
    using VirtualizedConstLocalFunctionWrapper =  VirtualizedConstLocalFunctionImpl< GridPart, Value, true /* wrapper */ >;

    //! VirtualizedConstLocalFunction stores the ConstLocalFunction as a object
    template< class GridPart, class Value >
    using VirtualizedConstLocalFunction =  VirtualizedConstLocalFunctionImpl< GridPart, Value, false /* wrapper */ >;

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_FUNCTION_VIRTUALIZEDGRIDFUNCTION_HH
