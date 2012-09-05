#ifndef DUNE_FEM_SHAPEFUNCTIONSET_LOCALFUNCTIONS_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_LOCALFUNCTIONS_HH

#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/shapefunctionset.hh>

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
      typedef typename Traits::Domain DomainType;
      typedef typename DomainType::value_type DomainFieldType;
      static const int dimDomain = DomainType::dimension;
      
      typedef typename Traits::Range RangeType;
      typedef typename RangeType::value_type RangeFieldType;
      static const int dimRange = RangeType::dimension;

      typedef FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange > FunctionSpaceType;
    };



    // LocalFunctionsShapeFunctionSet
    // ------------------------------

    template< class LocalBasis >
    class LocalFunctionsShapeFunctionSet
    : public ShapeFunctionSet< typename LocalFunctionsShapeFunctionSetTraits< LocalBasis >::FunctionSpaceType,
                               LocalFunctionsShapeFunctionSet< LocalBasis >
                             >
    {
      // this type
      typedef LocalFunctionsShapeFunctionSet< LocalBasis > ThisType;
      // traits class
      typedef LocalFunctionsShapeFunctionSetTraits< LocalBasis > Traits;
      // base type
      typedef ShapeFunctionSet< Traits, ThisType > BaseType;

    public:
      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
      typedef typename BaseType::RangeType RangeType;
      typedef typename BaseType::JacobianRangeType JacobianRangeType;

      explicit LocalFunctionsShapeFunctionSet ( const GeometryType &type, 
                                                const LocalBasis &localBasis )
      : type_( type ),
        localBasis_( localBasis )
      {
        values_.reserve( size() );
        jacobians_.reserve( size() );
      }

      GeometryType type () const { return type_; }

      std::size_t size () const { return localBasis_.size(); }

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

      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor &f ) const
      {
        DUNE_THROW( NotImplemented, "Method hessianEach not implemented" );
      }

    private:
      template< class T, class Functor >
      static void callFunctor ( const std::vector< T > &v, Functor &f )
      {
        typedef typename std::vector< T >::const_iterator Iterator;
        std::size_t i = 0;
        for( Iterator it = v.begin(); it != v.end(); ++it )
          f( i++, *it );
      }

      GeometryType type_;
      LocalBasis localBasis_;
      mutable std::vector< RangeType > values_;
      mutable std::vector< JacobianRangeType > jacobians_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_LOCALFUNCTIONS_HH
