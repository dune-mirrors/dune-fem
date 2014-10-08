#ifndef DUNE_FEM_SPACE_RANNACHERTUREK_LOCALINTERPOLATION_HH
#define DUNE_FEM_SPACE_RANNACHERTUREK_LOCALINTERPOLATION_HH

#include <cstddef>

#include <utility>
#include <vector>

#include <dune/common/fvector.hh>


namespace Dune
{

  namespace Fem
  {

    // RannacherTurekLocalInterpolation
    // --------------------------------

    template< class BasisFunctionSet, class LocalInterpolation >
    class RannacherTurekLocalInterpolation
    {
      typedef RannacherTurekLocalInterpolation< BasisFunctionSet, LocalInterpolation > ThisType;

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
      explicit RannacherTurekLocalInterpolation ( const BasisFunctionSetType &basisFunctionSet,
                                                  const LocalInterpolationType &localInterpolation = LocalInterpolationType() )
        : basisFunctionSet_( basisFunctionSet ),
          localInterpolation_( localInterpolation )
      {}

      RannacherTurekLocalInterpolation ( const ThisType & ) = default;

      RannacherTurekLocalInterpolation ( ThisType &&other )
        : basisFunctionSet_( std::move( other.basisFunctionSet_ ) ),
          localInterpolation_( std::move( other.localInterpolation_ ) )
      {}

      RannacherTurekLocalInterpolation &operator= ( const ThisType & ) = default;

      RannacherTurekLocalInterpolation &operator= ( ThisType &&other )
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
      void apply ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const;

    protected:
      const LocalInterpolationType &localInterpolation () const { return localInterpolation_; }

    private:
      BasisFunctionSet basisFunctionSet_;
      LocalInterpolationType localInterpolation_;
    };



    // Implementation of RannacherTurekLocalInterpolation::LocalFunctionWrapper
    // ------------------------------------------------------------------------

    template< class LocalInterpolation, class RangeVector >
    template< class LocalFunction >
    struct RannacherTurekLocalInterpolation< LocalInterpolation, RangeVector >::LocalFunctionWrapper
    {
      LocalFunctionWrapper ( const LocalFunction &localFunction, size_type component )
      : localFunction_( localFunction ),
        component_( component )
      {}

      template< class RangeFieldType >
      void evaluate ( const typename LocalFunction::DomainType &x, Dune::FieldVector< RangeFieldType, 1 > &y ) const
      {
        typename LocalFunction::RangeType z;
        localFunction().evaluate( x, z );
        y = z[ component() ];
      }

    private:
      const LocalFunction &localFunction () const { return localFunction_; }

      size_type component () const { return component_; }

      const LocalFunction &localFunction_;
      size_type component_;
    };



    // Implementation of RannacherTurekLocalInterpolation
    // --------------------------------------------------


    template< class LocalInterpolation, class RangeVector >
    template< class LocalFunction, class LocalDofVector >
    inline void RannacherTurekLocalInterpolation< LocalInterpolation, RangeVector >
      ::apply ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
    {
      std::vector< RangeFieldType > phi;

      for( int k = 0; k < dimRange; ++k )
      {
        LocalFunctionWrapper< LocalFunction > localFunctionWrapper( localFunction, k );
        localInterpolation().interpolate( localFunctionWrapper, phi );

        const size_type size = phi.size();
        for( size_type i = 0; i < size; ++i )
          localDofVector[ i*dimRange + k ] = phi[ i ];
      }
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_RANNACHERTUREK_LOCALINTERPOLATION_HH
