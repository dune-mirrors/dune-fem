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
        struct Traits
        {
          typedef typename LocalFunction::DomainType DomainType;
          typedef typename LocalFunction::RangeType RangeType;
        };
        typedef typename LocalFunction::DomainType DomainType;
        typedef typename LocalFunction::RangeType RangeType;

        LocalFunctionWrapper ( const LocalFunction &lf, const BasisFunctionSet &bset ) : lf_( lf ) {}

        template< class Arg >
        void evaluate ( const Arg &x, typename Traits::RangeType &y ) const
        {
          lf_.evaluate( x, y );
        }

      private:
        const LocalFunction &lf_;
      };


      // LocalFunctionWrapper
      // --------------------

      template< class LocalFunction, class Entity, class ShapeFunctionSet, class Transformation >
      struct LocalFunctionWrapper< LocalFunction, TransformedBasisFunctionSet< Entity, ShapeFunctionSet, Transformation > >
      {
        typedef TransformedBasisFunctionSet< Entity, ShapeFunctionSet, Transformation > BasisFunctionSetType;

        struct Traits
        {
          typedef typename LocalFunction::DomainType DomainType;
          typedef typename LocalFunction::RangeType RangeType;
        };
        typedef typename LocalFunction::DomainType DomainType;
        typedef typename LocalFunction::RangeType RangeType;

        LocalFunctionWrapper ( const LocalFunction &lf, const BasisFunctionSetType &bset ) : lf_( lf ), bset_( bset ) {}

        template< class Arg >
        void evaluate ( const Arg &x, typename Traits::RangeType &y ) const
        {
          typename Traits::RangeType help;
          lf_.evaluate( x, help );
          typename Transformation::InverseTransformationType transf( lf_.entity().geometry(), x );
          y = transf.apply( help );
        }

      private:
        const LocalFunction &lf_;
        const BasisFunctionSetType &bset_;
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

      template< class LocalFunction >
      using LocalFunctionWrapper = Impl::LocalFunctionWrapper< LocalFunction, BasisFunctionSetType >;

    public:
      explicit LocalFiniteElementInterpolation ( const BasisFunctionSetType &basisFunctionSet,
                                                 const LocalInterpolationType &localInterpolation = LocalInterpolationType() )
        : basisFunctionSet_( basisFunctionSet ),
        localInterpolation_( localInterpolation )
      {}

      template< class LocalFunction, class Dof >
      void operator() ( const LocalFunction &localFunction, std::vector< Dof > &localDofVector ) const
      {
        LocalFunctionWrapper< LocalFunction > wrapper( localFunction, basisFunctionSet() );
        localInterpolation().interpolate( wrapper, localDofVector );
      }

      template< class LocalFunction, class DiscreteFunction, template< class > class Assembly >
      void operator() ( const LocalFunction &localFunction, LocalContribution< DiscreteFunction, Assembly > &localContribution ) const
      {
        LocalFunctionWrapper< LocalFunction > wrapper( localFunction, basisFunctionSet() );
        localInterpolation().interpolate( wrapper, localContribution.localDofVector() );
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
