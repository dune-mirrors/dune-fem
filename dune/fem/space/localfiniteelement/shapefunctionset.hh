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

    template< class LocalBasis, int pSetId = -1 >
    class LocalFunctionsShapeFunctionSet
    {
      // this type
      typedef LocalFunctionsShapeFunctionSet< LocalBasis, pSetId > ThisType;
      // traits class
      typedef LocalFunctionsShapeFunctionSetTraits< LocalBasis > Traits;

    public:
      typedef typename Traits::FunctionSpaceType FunctionSpaceType;
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      static const int pointSetId = pSetId;

      explicit LocalFunctionsShapeFunctionSet ( const LocalBasis &localBasis )
        : localBasis_( localBasis )
      {
        values_.reserve( size() );
        jacobians_.reserve( size() );
        hessians_.resize( DomainType::dimension*(DomainType::dimension+1)/2. );
        for (unsigned int i=0;i<hessians_.size();++i)
          hessians_[i].reserve( size() );
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
        std::array<unsigned int, DomainType::dimension> multiIndex;
        std::fill(multiIndex.begin(),multiIndex.end(),0);
        unsigned int k = 0;
        for (unsigned int i=0;i<DomainType::dimension;++i)
        {
          multiIndex[i] = 1;
          for (unsigned int j=i;j<DomainType::dimension;++j)
          {
            multiIndex[j] += 1;
            localBasis_.partial(multiIndex,coordinate(x),hessians_[k]);
            multiIndex[j] -= 1;
            ++k;
          }
          multiIndex[i] -= 1;
        }
        callFunctor( hessians_, f );
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
      template< class T, class Functor >
      static void callFunctor ( const std::vector< std::vector<T> > &v, Functor f )
      {
        HessianRangeType h;
        for (unsigned int b=0;b<v[0].size();++b)
        {
          unsigned int k = 0;
          for (unsigned int i=0;i<DomainType::dimension;++i)
          {
            for (unsigned int j=i;j<DomainType::dimension;++j)
            {
              for (unsigned int r=0;r<RangeType::dimension;++r)
              {
                h[r][i][j] = v[b][k][r];
                h[r][j][i] = v[b][k][r];
              }
              ++k;
            }
          }
          f( b, h );
        }
      }

      const LocalBasis& localBasis_;
      mutable std::vector< RangeType > values_;
      mutable std::vector< JacobianRangeType > jacobians_;
      mutable std::vector< std::vector< RangeType > > hessians_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_LOCALFUNCTIONS_HH
