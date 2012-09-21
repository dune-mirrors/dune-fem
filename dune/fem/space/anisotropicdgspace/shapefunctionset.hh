#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SHAPEFUNCTIONSET_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SHAPEFUNCTIONSET_HH

// C++ includes
#include <cassert>

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/forloop.hh>

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>

// dune-fem includes
#include <dune/fem/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/legendre.hh>
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
    typedef LegendreShapeFunctionSetProvider< FunctionSpace, int maxOrder > ThisType;
    static const int storageSize = maxOrder+1;

  public:
    typedef FunctionSpace FunctionSpaceType;
    typedef LegendreShapeFunctionSet< FunctionSpaceType > LegendreShapeFunctionSetType;

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
      Dune::ForLoop< Initialize, 0, storageSize >::apply( shapeFunctionSets_);
    }

  public:
    ~LegendreShapeFunctionSetProvider ()
    {
      for( int i = 0; i < storageSize; ++i )
      {
        if( shapeFunctionSets_[ order ] )
          delete shapeFunctionSets_[ order ];
        shapeFunctionSets_[ order ] = nullptr;
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

    Dune::array< const LegendreShapeFunctionSet *, storageSize > shapeFunctionSets_;
  };



  // DefaultShapeFunctionSet
  // -----------------------

  /*
   * A class for constructing a shape function set.
   */
  template< class FunctionSpace, int maxOrder >
  struct DefaultShapeFunctionSet
  {
    // export template parameter
    typedef FunctionSpace FunctionSpaceType;

  private:
    typedef typename FunctionSpaceType::RangeType RangeType;
    static const int dimension = FunctionSpaceType::dimDomain;
    typedef typename Dune::ToScalarFunctionSpace< FunctionSpaceType >::Type ScalarFunctionSpace;
    typedef Dune::Fem::LegendreShapeFunctionSet< 
        typename ToLocalFunctionSpace< ScalarFunctionSpace, 1 >::Type
      > LegendreShapeFunctionSetType;
    typedef typename MakeTuple< const LegendreShapeFunctionSetType &, dimension >::Type LegendreShapeFunctionSetTupleType;
    typedef Dune::Fem::TensorProductShapeFunctionSet< ScalarFunctionSpace, LegendreShapeFunctionSetTupleType > ScalarShapeFunctionSet;

  public:
    // multi index type
    typedef MultiIndexSet< dimension, maxOrder >::MultiIndexType MultiIndexType;
    // implementation type
    typedef VectorialShapeFunctionSet< ScalarShapeFunctionSet, RangeType > ImplementationType;

    // return shape function set
    static ImplementationType create ( const MultiIndexType &multiIndex )
    {
      LegendreShapeFunctionSetTupleType tuple( multiIndex );
      ScalarShapeFunctionSet scalarShapeFunctionSet( tuple );
      return ImplementationType( scalarShapeFunctionSet );
    }
  };



  // ShapeFunctionSet
  // ----------------

  template< class FunctionSpace, int maxOrder, class Implementation = DefaultShapeFunctionSet< FunctionSpace, maxOrder > >
  class ShapeFunctionSet
  : public Dune::Fem::ShapeFunctionSet< FunctionSpace, ShapeFunctionSet< FunctionSpace, maxOrder, Implementation > >
  {
    typedef ShapeFunctionSet< FunctionSpace, maxOrder, Implementation > ThisType;
    typedef Dune::Fem::ShapeFunctionSet< FunctionSpace, ThisType > BaseType;

  public:
    typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
    typedef typename BaseType::RangeType RangeType;
    typedef typename BaseType::JacobianRangeType JacobianRangeType;
    typedef typename BaseType::HessianRangeType HessianRangeType;

  private:
    static const int dimension = FunctionSpaceType::dimDomain;

    typedef Implementation< FunctionSpaceType, maxOrder > ImplementationType;

  public:
    typedef typename ImplementationType::MultiIndexType MultiIndexType;

    ShapeFunctionSet ( const Dune::GeometryType &type, const MultiIndexType &multiIndex )
    : implementation_( multiIndex )
    {
      assert( type == ThisType::type() );
    }

    static Dune::GeometryType type ()
    {
      return implementation().type();
    }

    std::size_t size () const
    {
      return implementation().type();
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
    static ImplementationType ( const MultiIndexType &multiIndex ) const;

    ImplementationType implementation_;
  };

} // namespace AnisotropicDG 

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SHAPEFUNCTIONSET_HH
