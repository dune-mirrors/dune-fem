#ifndef DUNE_FEM_SPACE_FOURIER_BASISFUNCTIONSET_HH
#define DUNE_FEM_SPACE_FOURIER_BASISFUNCTIONSET_HH

#include <cassert>
#include <cstddef>

#include <dune/common/static_assert.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/fem/quadrature/quadrature.hh>
#include <dune/fem/space/basisfunctionset/functor.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/version.hh>

#include "basisfunctions.hh"


namespace Dune
{

  namespace Fem
  {

    // FourierBasisFunctionSet
    // -----------------------

    template< class Entity, class Range, int order >
    class FourierBasisFunctionSet
    {
      typedef FourierBasisFunctionSet< Entity, Range, order > ThisType;

      dune_static_assert( Range::dimension == 1, "FunctionSpace must be scalar (i.e., dimRange = 1)." );

    public:
      typedef Entity EntityType;

      typedef typename EntityType::Geometry GeometryType;

    private:
      typedef typename GeometryType::ctype DomainFieldType;
      static const int dimDomain = GeometryType::coorddimension;
      
    public:
      typedef std::size_t SizeType;

      typedef FunctionSpace< DomainFieldType, typename Range::value_type, dimDomain, Range::dimension > FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef Dune::ReferenceElement< typename DomainType::value_type, DomainType::dimension > ReferenceElementType;

    protected:
      typedef FourierBasisFunctions< FunctionSpaceType, order > BasisFunctionsType;

    public:
      /////////////////////////////////
      // Construction, copying, etc. //
      /////////////////////////////////

      FourierBasisFunctionSet ()
      : entity_( nullptr ), 
        geometry_( nullptr ),
        basisFunctions_()
      {}

      FourierBasisFunctionSet ( const ThisType &other )
      {
        *this = other;
      }

      const ThisType &operator= ( const ThisType &other )
      {
        entity_ = other.entity_;
        if( geometry_ )
          delete geometry_;
        if( entity_ )
          geometry_ = new GeometryType( entity_->geometry() );
      }

      ~FourierBasisFunctionSet ()
      {
        if( geometry_ )
          delete geometry_;
        geometry_ = nullptr;
      }


      //////////////////////////////////////////
      // Basis function set interface methods //
      //////////////////////////////////////////

      SizeType size () const { return basisFunctions().size(); }

      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
      {
        DomainType y = geometry().global( coordinate( x ) );
        FunctionalAxpyFunctor< RangeType, DofVector > f( valueFactor, dofs );
        basisFunctions().evaluateEach( y, f );
      }

      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        DomainType y = geometry().global( coordinate( x ) );
        FunctionalAxpyFunctor< JacobianRangeType, DofVector > f( jacobianFactor, dofs );
        basisFunctions().jacobianEach( y, f );
      }

      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor,
                  const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const
      {
        axpy( x, valueFactor, dofs );
        axpy( x, jacobianFactor, dofs );
      }

      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        const typename GeometryType::LocalCoordinate &y = coordinate( x );
        DomainType z = geometry().global( y );
        AxpyFunctor< DofVector, RangeType > f( dofs, value );
        basisFunctions().evaluateEach( z, f );
      }

      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        DomainType y = geometry().global( coordinate( x ) );
        AssignFunctor< RangeArray > f( values );
        basisFunctions().evaluateEach( y, f );
      }

      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        DomainType y = geometry().global( coordinate( x ) );
        AxpyFunctor< DofVector, JacobianRangeType > f( dofs, jacobian );
        basisFunctions().jacobianEach( y, f );
      }

      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        DomainType y = geometry().global( coordinate( x ) );
        AssignFunctor< JacobianRangeArray > f( jacobians );
        basisFunctions().jacobianEach( y, f );
      }

      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const
      {
        DomainType y = geometry().global( coordinate( x ) );
        AxpyFunctor< DofVector, HessianRangeType > f( dofs, hessian );
        basisFunctions().hessianEach( y, f );
      }

      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        DomainType y = geometry().global( coordinate( x ) );
        AssignFunctor< HessianRangeArray > f( hessians );
        basisFunctions().hessiansEach( y, f );
      }

      const EntityType &entity () const
      {
        assert( entity_ );
        return *entity_;
      }

      const GeometryType &geometry () const
      {
        assert( geometry_ );
        return *geometry_;
      }

      DUNE_VERSION_DEPRECATED(1,4,remove)
      Dune::GeometryType type () const { return geometry().type(); }

      const ReferenceElementType &referenceElement () const
      {
        return Dune::ReferenceElements< DomainFieldType, dimDomain >::general( geometry().type() );
      }


      ///////////////////////////
      // Non-interface methods //
      ///////////////////////////

      void initialize ( const Entity &entity ) const
      {
        entity_ = &entity;
        if( geometry_ )
          delete geometry_;
        geometry_ = new GeometryType( entity_->geometry() );
      }

    protected:
      const BasisFunctionsType basisFunctions () const { return basisFunctions_; }

    private:
      mutable const EntityType *entity_;
      mutable const GeometryType *geometry_;
      BasisFunctionsType basisFunctions_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FOURIER_BASISFUNCTIONSET_HH
