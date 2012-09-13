#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_LEGENDRESHAPEFUNCTIONSETS_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_LEGENDRESHAPEFUNCTIONSETS_HH

// C++ includes
#include <cassert>
#include <cstddef>

// dune-common includes
#include <dune/common/array.hh>
#include <dune/common/forloop.hh>
#include <dune/common/fvector.hh>
#include <dune/common/nullptr.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/dgspace/legendredgbasefunctions.hh>
#include <dune/fem/space/shapefunctionset/shapefunctionset.hh>

/**
  @file
  @author Christoph Gersbacher
  @brief Please doc me
*/


namespace AnisotropicDG
{

  // NumScalarShapeFunctions
  // -----------------------

  /**
   * \brief Provides dynamic access to the number of scalar 
   *        Legendre shape functions up to a given order
   *
   *  \tparam  dimension  Spatial dimension
   *  \tparam  maxOrder   Maximal polynomial order requested
   */
  template < int dimension, int maxOrder >
  class NumScalarShapeFunctions
  {
    typedef NumScalarShapeFunctions< dimension, maxOrder > ThisType;

    template< int order >
    struct Initialize
    {
      static void apply ( Dune::array< int, maxOrder+1 > &sizes )
      {
        dune_static_assert( order >= 0 && order <= maxOrder, "Invalid template parameter" );
        sizes[ order ] = Dune::Fem::NumLegendreBaseFunctions< order, dimension >::numBaseFct;
      }
    };

  protected:
    NumScalarShapeFunctions ()
    {
      Dune::ForLoop< Initialize, 0, maxOrder+1 >::apply( sizes_ );
    }

    static ThisType &instance ()
    {
      static ThisType instance_;
      return instance_;
    }
  
  public:
    /* \brief return number of scalar Legendre shape functions
     *
     * \note Legendre shape functions are only available on cubic reference elements.
     * \note type.dim() must coincide with template parameter dimension.
     */
    static std::size_t size ( const Dune::GeometryType &type, const int order ) 
    {
      assert( type.isCube() );
      assert( type.dim() == dimension );
      assert( order >= 0 && order <= maxOrder );
      return instance().sizes_[ order ];
    }

  private:
    NumScalarhapeFunctions ( const ThisType &other );
    ThisType &operator= ( const ThisType &other );

    Dune::array< std::size_t, maxOrder+1 > sizes_;
  };



  // SimpleShapeFunctionSet
  // ----------------------

  /**
   * \brief A very simple shape function set
   *
   * \tparam  Storage  An implementation of Dune::Fem::StorageInterface
   */
  template< class Storage >
  struct SimpleShapeFunctionSet
  {
    typedef typename StorageType::DomainType DomainType;
    typedef typename StorageType::RangeType RangeType;

    SimpleShapeFunctionSet ( const StorageType &storage )
    : storage_( &storage )
    { }

    template< int diffOrd >
    void evaluate( int i, const Dune::FieldVector< int, diffOrd > &diffVar,
                   const DomainType &x, RangeType &phi ) const
    {
      storage_->evaluate( i, diffVar, x, phi );
    }

  private:
    const StorageType *storage_;
  };



  // SimpleShapeFunctionsProvider
  // ----------------------------

  /**
   * \brief A class providing dynamic access to instances of SimpleShapeFunctionSet
   *        up to a given order
   *
   * \tparam  FunctionSpace  Scalar function space
   * \tparam  maxOrder       Maximal polynomial order requested
   * \tparam  Storage        A type of Dune::Fem::StorageInterface
   */
  template< class FunctionSpace, int maxOrder, template< class > class Storage = SimpleStorage >
  class SimpleShapeFunctionSetProvider
  {
    typedef SimpleShapeFunctionSetProvider< FunctionSpace, maxOrder, Storage > ThisType;

  public:
    typedef Storage< ScalarFunctionSpaceType > StorageType;

  private:
    template< int order >
    struct Initialize
    {
      // add scalar shape function set to storage
      static void apply ( const Dune::GeometryType &geometryType, array< StorageType, maxOrder+1 > &storage )
      {
        typedef LegendreDGBaseFunctionFactory< ScalarFunctionSpace, order > ScalarFactoryType;
        ScalarFactoryType factory( geometryType );
        storage[ order ] = StorageType( factory );
      }
    };

  protected:
    friend class DefaultSingletonFactory< const Dune::GeometryType, ThisType >;

    SimpleShapeFunctionSetProvider ( const Dune::GeometryType &type )
    {
      ForLoop< Initialize, 0, maxOrder+1 >::apply( type, storage_ );
    }

    static inline ThisType &instance ( const Dune::GeometryType &type )
    {
      typedef SingletonList< const Dune::GeometryType, ThisType > SingletonProviderType;
      return SingletonProviderType::getObject( type );
    }

  public:
    static SimpleShapeFunctionSetType ( const Dune::GeometryType &type, const int order )
    {
      assert( type.isCube() );
      assert( type.dim() == dimension );
      assert( order >= 0 && order <= maxOrder );
      return SimpleShapeFunctionSetType( instance( type ).storage_[ order ] );
    }

  private:
    SimpleShapeFunctionSetProvider ( const ThisType & );
    ThisType &operator= ( const ThisType & );

    array< StorageType, maxOrder+1 > storage_;
  };



  // ScalarShapeFunctionSetEvaluator
  // -------------------------------

  template< class FunctionSpace, int maxOrder >
  class ScalarShapeFunctionSetEvaluator
  {
    typedef ScalarShapeFunctionSetEvaluator< Functionspace, maxOrder > ThisType;

    typedef SimpleShapeFunctionsProvider< FunctionSpace, maxOrder > ScalarShapeFunctionsProviderType;
    typedef typename ScalarShapeFunctionsProviderType::StorageType StorageType;

  public:
    typedef typename FunctionSpace FunctionSpaceType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

    ScalarShapeFunctionSetEvaluator ( const Dune::GeometryType &type, const int order )
    : type_( type ),
      order_( order )
    { }

    std::size_t size () const
    {
      return NumScalarBaseFunctions< dimension, maxOrder >::size( type(), order() );
    }

    void evaluate( const int i, const DomainType &x, RangeType &value ) const
    {
      assert( i >= 0 && i < size() );
      Dune::FieldVector< int, 0 > diffVariable;
      shapeFunction( i ).evaluate( diffVariable, x, value );
    }

    void jacobian ( const int i, const DomainType &x, JacobianRangeType &jacobian ) const
    {
      assert( i >= 0 && i < size() );
      Dune::FieldVector< int, 1 > diffVariable;
      int &j = diffVariable[ 0 ];
      for( j = 0; j < dimension; ++j )
      {
        RangeType tmp;
        shapeFunction( i ).evaluate( diffVariable, x, tmp );
        jacobian[ 0 ][ j ] = tmp[ 0 ];
      }
    }

    void hessian ( const int baseFunction, const DomainType &x, HessianRangeType &jacobian ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method hessian() not implemented yet" );
    }

  protected:
    int order () const { return order_; }

    const Dune::GeometryType &type () const { return type_; }

    const StorageType &shapeFunction ( int i )
    {
      return ScalarShapeFunctionsProviderType::storage( type(), order() );
    }

  private:
    Dune::GeometryType type_;
    int order_;
  };



  // LegendreShapeFunctionSet
  // ------------------------

  template< class FunctionSpace, int maxOrder >
  class LegendreShapeFunctionSet
  : public Dune::Fem::ShapeFunctionSet< FunctionSpace, LegendreShapeFunctionSet< FunctionSpace, maxOrder > >
  {
    dune_static_assert( (FunctionSpace::dimRange == 1),
                        "FunctionSpace must be scalar (i.e., dimRange = 1)." );

    typedef LegendreShapeFunctionSet< FunctionSpace, maxOrder > ThisType;
    typedef Dune::Fem::ShapeFunctionSet< FunctionSpace, ThisType > BaseType;

    static const int dimension = FunctionSpace::dimDomain;
    typedef Dune::Fem::LegendreDGBaseFunction< FunctionSpace, maxOrder > ImplementationType;

  public:
    typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
    typedef typename BaseType::DomainType DomainType;
    typedef typename BaseType::RangeType RangeType;
    typedef typename BaseType::JacobianRangeType JacobianRangeType;
    typedef typename BaseType::HessianRangeType HessianRangeType;

    LegendreShapeFunctionSet ( const Dune::GeometryType &type, int p )
    : type_( type ),
      shapeFunctions_( numShapeFunctions, nullptr )
    {
      for( int i = 0; i < numShapeFunctions; ++i )
        shapeFunctions_[ i ] = new ImplementationType( i );
    }
    
    ~LegendreShapeFunctionSet ()
    {
      for( int i = 0; i < numShapeFunctions; ++i )
      {
        if( shapeFunctions_[ i ] )
          delete shapeFunctions_[ i ];
      }
    }

    Dune::GeometryType type () const
    {
      return type_;
    }

    std::size_t size () const
    {
      return NumScalarBaseFunctions< dimension, maxOrder >::size( type(), order() );
    }

    template< class Point, class Functor >
    void evaluateEach ( const Point &x, Functor functor ) const
    {
      for( std::size_t i = 0; i < size(); ++i )
      {
        RangeType val;
        evaluate( i, coordinate( x ), val );
        functor( i, val );
      }
    }

    template< class Point, class Functor >
    void jacobianEach ( const Point &x, Functor functor ) const
    {
      for( std::size_t i = 0; i < size(); ++i )
      {
        JacobianRangeType jac;
        jacobian( i, coordinate( x ), jac );
        functor( i, jac );
      }
    }

    template< class Point, class Functor >
    void hessianEach ( const Point &x, Functor functor ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method hessianEach() not implemented yet" );
    }

  private:
    Dune::GeometryType type_; 
    Dune::array< ImplementationType *, numShapeFunctions > shapeFunctions_;
  };

} // namespace AnisotropicDG 

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_LEGENDRESHAPEFUNCTIONSETS_HH
