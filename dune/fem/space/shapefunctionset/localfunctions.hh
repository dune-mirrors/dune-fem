#ifndef DUNE_FEM_SHAPEFUNCTIONSET_LOCALFUNCTIONS_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_LOCALFUNCTIONS_HH

#include <vector>

namespace Dune
{

  namespace Fem
  {

    // LocalFunctionsShapeFunctionSet
    // ------------------------------

    template< class LocalBasis >
    class LocalFunctionsShapeFunctionSet
    {
      typedef LocalFunctionsShapeFunctionSet< LocalBasis > ThisType;

    public:
      typedef LocalBasis LocalBasisType;

      typedef typename LocalBasis::Traits::Range RangeType;
      typedef typename LocalBasis::Traits::Jacobian JacobianRangeType;

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor &f ) const
      {
        localBasis_.evaluateFunction( coordinate( x ), values_ );
        assert( jacobians_.size() == size() );
        callFunctor( values_, f );
      }

      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor &f ) const
      {
        localBasis_.jacobianFunction( coordinate( x ), jacobians_ );
        assert( jacobians_.size() == size() );
        callFunctor( jacobians_, f );
      }

      std::size_t size () const { return localBasis_.size(); }

    private:
      template< class T, class Functor >
      static void callFunctor ( const std::vector< T > &v, Functor &f )
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
