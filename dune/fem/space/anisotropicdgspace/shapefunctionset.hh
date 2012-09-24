#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SHAPEFUNCTIONSET_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SHAPEFUNCTIONSET_HH

// C++ includes
#include <cassert>
#include <cstddef>

// dune-common includes

#include <dune/common/array.hh>
#include <dune/common/densevector.hh>
#include <dune/common/forloop.hh>
#include <dune/common/static_assert.hh>

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>

// dune-fem includes
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/dgspace/legendredgbasefunctions.hh>
#include <dune/fem/space/shapefunctionset/legendre.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/tensorproduct.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

// local includes
#include "multiindexset.hh"
#include "utility.hh"

/**
  @file
  @author Christoph Gersbacher
  @brief Please doc me
*/


namespace AnisotropicDG 
{

  // NumShapeFunctions
  // -----------------

  /**
   * \brief Provide number of shape functions for AnisotropicDGSpace
   *
   * \tparam  dimension  grid dimension
   *
   * \tparam  maxOrder   maximum polynomal order
   */
  template < int dimension, int maxOrder >
  class NumShapeFunctions
  {
    // this type
    typedef NumShapeFunctions< dimension, maxOrder > ThisType;

    template< int order >
    struct Initialize
    {
      static void apply ( Dune::array< std::size_t, maxOrder+1 > &sizes )
      {
        dune_static_assert( order >= 0 && order <= maxOrder, "Invalid template parameter " );
        sizes[ order ] = Dune::Fem::NumLegendreBaseFunctions< order, dimension >::numBaseFct;
      }
    };

  protected:
    //! \brief constructor
    NumShapeFunctions ()
    {
      Dune::ForLoop< Initialize, 0, maxOrder >::apply( sizes_ );
    }

    //! \brief get singleton
    static ThisType &instance ()
    {
      static ThisType instance_;
      return instance_;
    }

  public:
    //! \brief return max number of shape functions
    static std::size_t max ()
    {
      std::size_t min = 1;
      for( int i = 0; i < dimension; ++i )
        min *= instance().sizes_[ maxOrder ];
      return min ;
    }

    //! \brief return min number of shape functions
    static std::size_t min ()
    {
      std::size_t min = 1;
      for( int i = 0; i < dimension; ++i )
        min *= instance().sizes_[ 0 ];
      return min ;
    }

    //! \brief return number of shape functions for given multi index
    template< class Implementation >
    static std::size_t count ( const Dune::DenseVector< Implementation > &multiIndex )
    {
      assert( multiIndex.size() == dimension );
      std::size_t count = 1;
      for( int i = 0; i < dimension; ++i )
        count *= instance().sizes_[ multiIndex[ i ] ];
      return count;
    }

  private:
    //forbid copy constructor
    NumShapeFunctions ( const ThisType &other );
    // forbid assignment operator
    ThisType &operator= ( const ThisType &other );

    Dune::array< std::size_t, maxOrder+1 > sizes_;
  };



  // LegendreShapeFunctionSetProvider
  // --------------------------------

  /*
   * \brief Singleton class providing access to instances of LegendreShapeFunctionSet
   *
   * There is only one public method on this class:
\code
    static const LegendreShapeFunctionSetType &get ( const int order );
\endcode
    *
    * Example usage:
\code
    int order;
    ...
    typedef typename LegendreShapeFunctionSetProviderType::LegendreShapeFunctionSetType 
      LegendreShapeFunctionSetType;
    const LegendreShapeFunctionSetType &legendreShapeFunctionSet 
      = LegendreShapeFunctionSetProviderType::get( order );
\endcode
  */
  template< class FunctionSpace, int maxOrder >
  class LegendreShapeFunctionSetProvider
  {
    typedef LegendreShapeFunctionSetProvider< FunctionSpace, maxOrder > ThisType;
    static const int storageSize = maxOrder+1;

  public:
    typedef FunctionSpace FunctionSpaceType;
    typedef Dune::Fem::LegendreShapeFunctionSet< FunctionSpaceType > LegendreShapeFunctionSetType;

  private:
    template< int order >
    struct Initialize
    {
      static void apply ( const Dune::array< LegendreShapeFunctionSetType *, storageSize > &shapeFunctions )
      {
        dune_static_assert( order >= 0 && order < storageSize, "Invalid template parameter" );
        shapeFunctions[ order ] = new LegendreShapeFunctionSetType( order );
      }
    };

  protected:
    LegendreShapeFunctionSetProvider ()
    {
      Dune::ForLoop< Initialize, 0, maxOrder >::apply( shapeFunctionSets_);
    }

  public:
    ~LegendreShapeFunctionSetProvider ()
    {
      for( int i = 0; i < storageSize; ++i )
      {
        if( shapeFunctionSets_[ i ] )
          delete shapeFunctionSets_[ i ];
        shapeFunctionSets_[ i ] = nullptr;
      }
    }

    static const LegendreShapeFunctionSetType &get ( const int order )
    {
      return *( instance().shapeFunctionSets_[ order ] );
    }

  protected:
    static const ThisType &instance ()
    {
      static ThisType instance_;
      return instance_;
    }

  private:
    LegendreShapeFunctionSetProvider ( const ThisType & );
    ThisType &operator= ( const ThisType & );

    Dune::array< const LegendreShapeFunctionSetType *, storageSize > shapeFunctionSets_;
  };



  // ShapeFunctionSetTupleProvider
  // -----------------------------

  /*
   * \brief A class providing tuples of shape function sets
   *
   * A class providing tuples of one-dimensional, scalar shape function sets
   * to be used with Dune::Fem::TensorShapeFunctionSet
   *
   * Example usage:
\code
    const int maxOrder = MAXORDER;

    // create multi index
    typedef typedef ShapeFunctionSetTupleProvider< FunctionSpaceType, maxOrder >::MultiIndexType 
      MultiIndexType;
    MultiIndexType multiIndex( maxOrder );

    // create shape function set tuple
    typedef typedef ShapeFunctionSetTupleProvider< FunctionSpaceType, maxOrder >::ShapeFunctionSetTupleType
      ShapeFunctionSetTupleType;
    const ShapeFunctionSetTupleType shapeFunctionSetTuple
      = ShapeFunctionSetTupleProvider< FunctionSpaceType, maxOrder >::create( multiIndex );
\endcode
   */
  template< class FunctionSpace, int maxOrder >
  struct ShapeFunctionSetTupleProvider
  {
    // export template parameter function space
    typedef FunctionSpace FunctionSpaceType;

  private:
    // dimension
    static const int dimension = FunctionSpaceType::dimDomain;

    // scalar function space
    typedef typename Dune::Fem::ToScalarFunctionSpace< FunctionSpaceType >::Type ScalarFunctionSpace;

    // provider for simple shape function sets 
    typedef LegendreShapeFunctionSetProvider< 
        typename Dune::Fem::ToLocalFunctionSpace< ScalarFunctionSpace, 1 >::Type, maxOrder
      > ShapeFunctionSetProviderType;

    // simple shape function sets
    typedef typename ShapeFunctionSetProviderType::LegendreShapeFunctionSetType ShapeFunctionSetType;

    // we need to apply a proxy here
    typedef Dune::Fem::ShapeFunctionSetProxy< ShapeFunctionSetType > ShapeFunctionSetProxyType;

  public:
    //! \brief type of shape function set tuple
    typedef typename MakeTuple< ShapeFunctionSetProxyType, dimension >::Type ShapeFunctionSetTupleType;

    //! \brief multi index type
    typedef typename MultiIndexSet< dimension, maxOrder >::MultiIndexType MultiIndexType;

  protected:
    template< int i >
    struct Create
    {
      static void apply ( ShapeFunctionSetTupleType &shapeFunctionSetTuple, const MultiIndexType &multiIndex )
      {
        assert( MultiIndexSetType::contains( multiIndex ) );
        Dune::get< i >( shapeFunctionSetTuple )
          = &( ShapeFunctionSetProviderType::get( multiIndex[ i ] ) );
      }
    };

  public:
    //! \brief create shape function set tuple from multi index
    static ShapeFunctionSetTupleType create ( const MultiIndexType &multiIndex )
    {
      ShapeFunctionSetTupleType shapeFunctionSetTuple;
      Dune::ForLoop< Create, 0, dimension >::apply( shapeFunctionSetTuple, multiIndex );
      return shapeFunctionSetTuple;
    }
  };



  // ScalarShapeFunctionSet
  // ----------------------

  /*
   * \brief Provices an anisotropic scalar Legendre shape function set
   */
  template< class FunctionSpace, int maxOrder >
  struct ScalarShapeFunctionSet
  {
    typedef FunctionSpace FunctionSpaceType;

  private:
    dune_static_assert( (FunctionSpaceType::dimRange == 1),
                        "FunctionSpace must be scalar (i.e., dimRange = 1)." );

    typedef ScalarShapeFunctionSet< FunctionSpaceType, maxOrder > ThisType;

  private:
    typedef ShapeFunctionSetTupleProvider< FunctionSpaceType, maxOrder > ShapeFunctionSetTupleProviderType;
    typedef typename ShapeFunctionSetTupleProviderType::ShapeFunctionSetTupleType ShapeFunctionSetTupleType;

  protected:
    typedef Dune::Fem::TensorProductShapeFunctionSet< 
        FunctionSpaceType, ShapeFunctionSetTupleType
      > ImplementationType;

  public:
    typedef typename ShapeFunctionSetTupleProviderType::MultiIndexType MultiIndexType;

    ScalarShapeFunctionSet ( const MultiIndexType &multiIndex )
    : implementation_( ShapeFunctionSetTupleType::create( multiIndex ) )
    {
      assert(( NumShapeFunctions< FunctionSpaceType::dimDomain, maxOrder >::count( multiIndex ) == size() )); 
    }

    Dune::GeometryType type () const { return implementation().type(); }

    std::size_t size () const
    {
      return implementation().size();
    }

    template< class Point, class Functor >
    void evaluateEach ( const Point &x, Functor functor ) const
    {
      return implementation().evaluateEach( x, functor );
    }

    template< class Point, class Functor >
    void jacobianEach ( const Point &x, Functor functor ) const
    {
      return implementation().jacobianEach( x, functor );
    }

    template< class Point, class Functor >
    void hessianEach ( const Point &x, Functor functor ) const
    {
      return implementation().hessianEach( x, functor );
    }

  protected:
    const ImplementationType &implementation () const
    {
      return implementation_;
    }

  private:
    ImplementationType implementation_;
  };



  // ShapeFunctionSet
  // ----------------

  template< class FunctionSpace, int maxOrder >
  class ShapeFunctionSet
  {
    typedef ShapeFunctionSet< FunctionSpace, maxOrder > ThisType;

  public:
    typedef FunctionSpace FunctionSpaceType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  private:
    typedef typename Dune::Fem::ToScalarFunctionSpace< FunctionSpaceType >::Type ScalarFunctionSpace;
    typedef ScalarShapeFunctionSet< ScalarFunctionSpace, maxOrder > ScalarShapeFunctionSetType;

  protected:
    typedef Dune::Fem::VectorialShapeFunctionSet< ScalarShapeFunctionSetType, RangeType > ImplementationType;

  public:
    typedef typename ScalarShapeFunctionSetType::MultiIndexType MultiIndexType;

    ShapeFunctionSet ( const Dune::GeometryType &type, const MultiIndexType &multiIndex )
    : implementation_( ScalarShapeFunctionSetType( multiIndex ) )
    {
      assert( type == ThisType::type() );
    }

    Dune::GeometryType type () { return implementation().type(); }

    std::size_t size () const { return implementation().size(); }

    template< class Point, class Functor >
    void evaluateEach ( const Point &x, Functor functor ) const
    {
      return implementation().evaluateEach( x, functor );
    }

    template< class Point, class Functor >
    void jacobianEach ( const Point &x, Functor functor ) const
    {
      return implementation().jacobianEach( x, functor );
    }

    template< class Point, class Functor >
    void hessianEach ( const Point &x, Functor functor ) const
    {
      return implementation().hessianEach( x, functor );
    }

  protected:
    const ImplementationType &implementation () const
    {
      return implementation_;
    }

  private:
    ImplementationType implementation_;
  };

} // namespace AnisotropicDG 

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SHAPEFUNCTIONSET_HH
