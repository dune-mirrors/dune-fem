#ifndef DUNE_FEM_SPACE_GENERICDISCRETE_LOCALINTERPOLATION_HH
#define DUNE_FEM_SPACE_GENERICDISCRETE_LOCALINTERPOLATION_HH

// C++ includes
#include <cstdlib>
#include <vector>

// dune-common includes
#include <dune/common/fvector.hh>


namespace Dune
{

  namespace Fem
  {

    // VectorialLocalInterpolation
    // ---------------------------

    template< class LocalInterpolation, class RangeVector >
    struct VectorialLocalInterpolation
    {
      typedef LocalInterpolation LocalInterpolationType;
      typedef RangeVector RangeType;

      typedef typename RangeType::field_type RangeFieldType;
      static const int dimRange = RangeType::dimension;

      typedef std::size_t size_type;

    private:
      template< class LocalFunction >
      struct LocalFunctionWrapper;

    public:
      explicit VectorialLocalInterpolation ( const LocalInterpolationType &localInterpolation = LocalInterpolationType() )
      : localInterpolation_( localInterpolation )
      {}
      
      template< class LocalFunction, class LocalDofVector >
      void operator() ( const LocalFunction &f, LocalDofVector &dofs ) const;

    protected:
      const LocalInterpolationType &localInterpolation () const { return localInterpolation_; }

    private:
      LocalInterpolationType localInterpolation_;
    };



    // Implementation of VectorialLocalInterpolation::LocalFunctionWrapper
    // -------------------------------------------------------------------

    template< class LocalInterpolation, class RangeVector >
    template< class LocalFunction >
    struct VectorialLocalInterpolation< LocalInterpolation, RangeVector >::LocalFunctionWrapper
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



    // Implementation of VectorialLocalInterpolation
    // ---------------------------------------------

     
    template< class LocalInterpolation, class RangeVector >
    template< class LocalFunction, class LocalDofVector >
    inline void VectorialLocalInterpolation< LocalInterpolation, RangeVector >
      ::operator() ( const LocalFunction &f, LocalDofVector &dofs ) const
    {
      std::vector< RangeFieldType > phi;
      
      for( int k = 0; k < dimRange; ++k )
      {
        LocalFunctionWrapper< LocalFunction > localFunctionWrapper( f, k );
        localInterpolation().interpolate( localFunctionWrapper, phi );
        
        const size_type size = phi.size();
        for( size_type i = 0; i < size; ++i )
        {
          dofs[ i*dimRange + k ] = phi[ i ];
        }
      }
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_GENERICDISCRETE_LOCALINTERPOLATION_HH
