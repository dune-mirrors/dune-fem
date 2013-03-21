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

    template< class Entity, class BasisFunctions >
    class FourierBasisFunctionSet
    {
      typedef FourierBasisFunctionSet< Entity, BasisFunctions > ThisType;

    public:
      typedef Entity EntityType;
      typedef BasisFunctions BasisFunctionsType;

      typedef typename EntityType::Geometry GeometryType;

      typedef std::size_t SizeType;

      typedef typename BasisFunctionsType::FunctionSpaceType FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef Dune::ReferenceElement< typename DomainType::value_type, DomainType::dimension > ReferenceElementType;

      FourierBasisFunctionSet ()
      : entity_( nullptr ),
        basisFunctions_( nullptr )
      {}

      FourierBasisFunctionSet ( const EntityType &entity, const BasisFunctionsType &basisFunctions )
      : entity_( &entity ),
        basisFunctions_( &basisFunctions )
      {}

      //////////////////////////////////////////
      // Basis function set interface methods //
      //////////////////////////////////////////

      int order () const { return basisFunctions().order(); }

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
        DomainType y = geometry().global( coordinate( x ) );
        FunctionalAxpyFunctor< RangeType, DofVector > fv( valueFactor, dofs );
        basisFunctions().evaluateEach( y, fv );
        FunctionalAxpyFunctor< JacobianRangeType, DofVector > fj( jacobianFactor, dofs );
        basisFunctions().jacobianEach( y, fj );
      }

      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        DomainType y = geometry().global( coordinate( x ) );
        AxpyFunctor< DofVector, RangeType > f( dofs, value );
        basisFunctions().evaluateEach( y, f );
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
        basisFunctions().hessianEach( y, f );
      }

      const EntityType &entity () const
      {
        assert( entity_ );
        return *entity_;
      }

      DUNE_VERSION_DEPRECATED(1,4,remove)
      Dune::GeometryType type () const { return geometry().type(); }

      const ReferenceElementType &referenceElement () const
      {
        return Dune::ReferenceElements< typename GeometryType::ctype, GeometryType::dimension >::general( geometry().type() );
      }

    protected:
      GeometryType geometry () const { return entity().geometry(); }

      const BasisFunctionsType basisFunctions () const
      {
        assert( basisFunctions_ );
        return *basisFunctions_;
      }

    private:
      const EntityType *entity_;
      const BasisFunctionsType *basisFunctions_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FOURIER_BASISFUNCTIONSET_HH
