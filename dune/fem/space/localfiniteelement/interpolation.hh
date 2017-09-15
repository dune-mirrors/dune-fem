#ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_INTERPOLATION_HH
#define DUNE_FEM_SPACE_LOCALFINITEELEMENT_INTERPOLATION_HH

#include <cstddef>

#include <utility>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/space/basisfunctionset/transformed.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Impl
    {

      // LocalFunctionWrapper
      // --------------------

      template< class LocalFunction, class BasisFunctionSet >
      struct LocalFunctionWrapper
      {
        typedef typename BasisFunctionSet::ShapeFunctionSetType ShapeFunctionSetType;

        struct Traits
        {
          typedef typename ShapeFunctionSetType::ScalarShapeFunctionSetType::DomainType DomainType;
          typedef typename ShapeFunctionSetType::ScalarShapeFunctionSetType::RangeType RangeType;
        };

        typedef typename Traits::DomainType DomainType;
        typedef typename Traits::RangeType RangeType;

        LocalFunctionWrapper ( const LocalFunction &lf, const BasisFunctionSet &bset, int component )
          : lf_( lf ), component_( component )
        {}

        void evaluate ( const DomainType &x, RangeType &y ) const
        {
          typedef MakeVectorialTraits< RangeType, typename ShapeFunctionSetType::RangeType > Traits;
          y = Traits::component( lf_.evaluate( x ), component_ );
        }

      private:
        const LocalFunction &lf_;
        int component_;
      };


      // LocalFunctionWrapper
      // --------------------

      template< class LocalFunction, class Entity, class ShapeFunctionSet, class Transformation >
      struct LocalFunctionWrapper< LocalFunction, TransformedBasisFunctionSet< Entity, ShapeFunctionSet, Transformation > >
      {
        typedef TransformedBasisFunctionSet< Entity, ShapeFunctionSet, Transformation > BasisFunctionSetType;

        struct Traits
        {
          typedef typename ShapeFunctionSet::ScalarShapeFunctionSetType::DomainType DomainType;
          typedef typename ShapeFunctionSet::ScalarShapeFunctionSetType::RangeType RangeType;
        };

        typedef typename Traits::DomainType DomainType;
        typedef typename Traits::RangeType RangeType;

        typedef Fem::MakeVectorialTraits< RangeType, typename ShapeFunctionSet::RangeType > MakeVectorialTraits;

        LocalFunctionWrapper ( const LocalFunction &lf, const BasisFunctionSetType &bset, typename MakeVectorialTraits::IndexType component )
          : lf_( lf ), bset_( bset ), component_( component )
        {}

        void evaluate ( const DomainType &x, RangeType &y ) const
        {
          typedef typename Transformation::InverseTransformationType Trafo;
          y = MakeVectorialTraits::component( Trafo( lf_.entity().geometry(), x ).apply( lf_.evaluate( x ) ), component_ );
        }

      private:
        const LocalFunction &lf_;
        const BasisFunctionSetType &bset_;
        typename MakeVectorialTraits::IndexType component_;
      };

    } // namespace Impl



    // LocalFiniteElementInterpolation
    // -------------------------------

    template< class BasisFunctionSet, class LocalInterpolation >
    class LocalFiniteElementInterpolation
    {
      typedef LocalFiniteElementInterpolation< BasisFunctionSet, LocalInterpolation > ThisType;

    public:
      typedef BasisFunctionSet BasisFunctionSetType;
      typedef LocalInterpolation LocalInterpolationType;

    private:
      typedef typename BasisFunctionSetType::FunctionSpaceType FunctionSpaceType;

      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
      static const int dimRange = FunctionSpaceType::dimRange;

      typedef std::size_t size_type;

      typedef typename BasisFunctionSetType::ShapeFunctionSetType ShapeFunctionSetType;
      typedef typename ShapeFunctionSetType::ScalarShapeFunctionSetType ScalarShapeFunctionSetType;

      typedef Fem::MakeVectorialTraits< typename ScalarShapeFunctionSetType::RangeType, typename ShapeFunctionSetType::RangeType > MakeVectorialTraits;

      template< class LocalFunction >
      using LocalFunctionWrapper = Impl::LocalFunctionWrapper< LocalFunction, BasisFunctionSetType >;

    public:
      explicit LocalFiniteElementInterpolation ( const BasisFunctionSetType &basisFunctionSet,
                                                 const LocalInterpolationType &localInterpolation = LocalInterpolationType() )
        : basisFunctionSet_( basisFunctionSet ),
          localInterpolation_( localInterpolation )
      {}

      template< class LocalFunction, class LocalDofVector >
      auto operator() ( const LocalFunction &localFunction, LocalDofVector &&localDofVector ) const
        -> void_t< decltype( localDofVector[ 0 ] = RangeFieldType( 0 ) ) >
      {
        std::vector< RangeFieldType > phi;

        const typename MakeVectorialTraits::IndexType end = MakeVectorialTraits::end();
        for( typename MakeVectorialTraits::IndexType k = MakeVectorialTraits::begin(); k != end; ++k )
        {
          LocalFunctionWrapper< LocalFunction > wrapper( localFunction, basisFunctionSet(), k );
          localInterpolation().interpolate( wrapper, phi );

          for( size_type i = 0, size = phi.size(); i < size; ++i )
            localDofVector[ i*MakeVectorialTraits::factor + MakeVectorialTraits::index( k ) ] = phi[ i ];
        }
      }

      BasisFunctionSetType basisFunctionSet () const { return basisFunctionSet_; }
      const LocalInterpolationType &localInterpolation () const { return localInterpolation_; }

    private:
      BasisFunctionSet basisFunctionSet_;
      LocalInterpolationType localInterpolation_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_INTERPOLATION_HH
