#ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_INTERPOLATION_HH
#define DUNE_FEM_SPACE_LOCALFINITEELEMENT_INTERPOLATION_HH

#include <cstddef>

#include <utility>
#include <vector>

#include <dune/common/fvector.hh>

namespace Dune
{

  namespace Fem
  {

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
      struct LocalFunctionWrapper;
      {
        struct Traits
        {
          typedef typename LocalFunction::RangeType RangeType;
        };

        LocalFunctionWrapper ( const LocalFunction &lf ) : lf_( lf ) {}

        template< class Arg >
        void evaluate ( const Arg &x, typename Traits::RangeType &y ) const
        {
          lf_.evaluate( x, y );
        }

      private:
        const LocalFunction &lf_;
      };

    public:
      explicit LocalFiniteElementInterpolation ( const BasisFunctionSetType &basisFunctionSet,
                                                 const LocalInterpolationType &localInterpolation = LocalInterpolationType() )
        : basisFunctionSet_( basisFunctionSet ),
          localInterpolation_( localInterpolation )
      {}

      template< class LocalFunction, class LocalDofVector >
      void operator() ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        LocalFunctionWrapper< LocalFunction > wrapper( localFunction );
        localInterpolation().interpolate( wrapper, localDofVector );
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
