#ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_LOCALINTERPOLATION_HH
#define DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_LOCALINTERPOLATION_HH

#include <cstddef>

#include <utility>
#include <vector>

#include <dune/common/fvector.hh>


namespace Dune
{

  namespace Fem
  {

    // BDMLocalInterpolation
    // ---------------------

    template< class BasisFunctionSet, class LocalInterpolation >
    class BDMLocalInterpolation
    {
      typedef BDMLocalInterpolation< BasisFunctionSet, LocalInterpolation > ThisType;

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

    public:
      explicit BDMLocalInterpolation ( const BasisFunctionSetType &basisFunctionSet,
                                       const LocalInterpolationType &localInterpolation = LocalInterpolationType() )
        : basisFunctionSet_( basisFunctionSet ),
          localInterpolation_( localInterpolation )
      {}

      BDMLocalInterpolation ( const ThisType & ) = default;

      BDMLocalInterpolation ( ThisType &&other )
        : basisFunctionSet_( std::move( other.basisFunctionSet_ ) ),
          localInterpolation_( std::move( other.localInterpolation_ ) )
      {}

      BDMLocalInterpolation &operator= ( const ThisType & ) = default;

      BDMLocalInterpolation &operator= ( ThisType &&other )
      {
        basisFunctionSet_ = std::move( other.basisFunctionSet_ );
        localInterpolation_ = std::move( other.localInterpolation_ );
        return *this;
      }

      BasisFunctionSetType basisFunctionSet () const
      {
        return basisFunctionSet_;
      }

      template< class LocalFunction, class LocalDofVector >
      void operator() ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        apply( localFunction, localDofVector );
      }

      template< class LocalFunction, class LocalDofVector >
      void apply ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        LocalFunctionWrapper< LocalFunction > wrapper( localFunction );
        localInterpolation().interpolate( wrapper, localDofVector );
      }

    protected:
      const LocalInterpolationType &localInterpolation () const { return localInterpolation_; }

    private:
      BasisFunctionSet basisFunctionSet_;
      LocalInterpolationType localInterpolation_;
    };


    // Implementation of BDMLocalInterpolation::LocalFunctionWrapper
    // -------------------------------------------------------------

    template< class BasisFunctionSet, class LocalInterpolation >
    template< class LocalFunction >
    struct BDMLocalInterpolation< BasisFunctionSet, LocalInterpolation >::LocalFunctionWrapper
    {
      struct Traits
      {
        typedef typename LocalFunction::RangeType RangeType;
      };

      LocalFunctionWrapper( const LocalFunction &lf ) : lf_( lf ) {}

      template< class Arg >
      void evaluate ( const Arg &x, typename LocalFunction::RangeType &y ) const
      {
        typename LocalFunction::RangeType help;
        lf_.evaluate( x, help );

        auto DF = lf_.entity().geometry().jacobianInverseTransposed( x );

        typename LocalFunction::RangeFieldType det = 1.0 / DF.determinant();

        y = 0;
        DF.usmtv( det, help, y );
      }

    private:
      const LocalFunction &lf_;

    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_LOCALINTERPOLATION_HH
