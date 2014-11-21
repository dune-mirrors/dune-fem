#ifndef DUNE_FEM_SHAPEFUNCTIONSET_LOCALFUNCTIONS_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_LOCALFUNCTIONS_HH

// C++ includes
#include <cstddef>
#include <vector>

// dune-common includes
#include <dune/common/exceptions.hh>

// dune-fem includes
#include <dune/fem/common/coordinate.hh>
#include <dune/fem/space/common/functionspace.hh>

namespace Dune
{

  namespace Fem
  {

    // LocalFunctionsShapeFunctionSetTraits
    // ------------------------------------

    template< class LocalBasis >
    class LocalFunctionsShapeFunctionSetTraits
    {

      typedef typename LocalBasis::Traits Traits;

    public:
      typedef typename Traits::DomainType DomainType;
      typedef typename DomainType::value_type DomainFieldType;
      static const int dimDomain = DomainType::dimension;

      typedef typename Traits::RangeType RangeType;
      typedef typename RangeType::value_type RangeFieldType;
      static const int dimRange = RangeType::dimension;

      typedef FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;
    };



    // LocalFunctionsShapeFunctionSet
    // ------------------------------

    template< class LocalBasis >
    class LocalFunctionsShapeFunctionSet
    {
      // this type
      typedef LocalFunctionsShapeFunctionSet< LocalBasis > ThisType;
      // traits class
      typedef LocalFunctionsShapeFunctionSetTraits< LocalBasis > Traits;

    public:
      typedef typename Traits::FunctionSpaceType FunctionSpaceType;
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      explicit LocalFunctionsShapeFunctionSet ( const LocalBasis &localBasis )
      : localBasis_( localBasis )
      {
        values_.reserve( size() );
        jacobians_.reserve( size() );
      }

      int order () const { return localBasis_.order(); }

      std::size_t size () const { return localBasis_.size(); }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor f ) const
      {
        localBasis_.evaluateFunction( coordinate( x ), values_ );
        assert( values_.size() == size() );
        callFunctor( values_, f );
      }

      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor f ) const
      {
        localBasis_.evaluateJacobian( coordinate( x ), jacobians_ );
        assert( jacobians_.size() == size() );
        callFunctor( jacobians_, f );
      }

      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor f ) const
      {
        DUNE_THROW( NotImplemented, "Method hessianEach not implemented" );
      }

    private:
      template< class T, class Functor >
      static void callFunctor ( const std::vector< T > &v, Functor f )
      {
        typedef typename std::vector< T >::const_iterator Iterator;
        std::size_t i = 0;
        for( Iterator it = v.begin(); it != v.end(); ++it )
          f( i++, *it );
      }

      LocalBasis localBasis_;
      mutable std::vector< RangeType > values_;
      mutable std::vector< JacobianRangeType > jacobians_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_LOCALFUNCTIONS_HH
