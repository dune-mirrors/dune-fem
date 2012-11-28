#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SHAPEFUNCTIONSETSTORAGE_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SHAPEFUNCTIONSETSTORAGE_HH

// C++ includes
#include <cassert>
#include <cstddef>
#include <vector>

// dune-common includes
#include <dune/common/nullptr.hh>
#include <dune/common/power.hh>

// dune-geometry
#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/storage/singletonlist.hh>

// local includes
#include "multiindexset.hh"
#include "shapefunctionset.hh"

/**
  @file
  @author Christoph Gersbacher
  @brief Please doc me.
*/


namespace AnisotropicDG
{

  // ShapeFunctionSetStorage
  // -----------------------

  /** 
   * \brief Storage class for shape function sets for anisotropic DG space
   *
   * \tparam  FunctionSpace  Scalar value function space
   * \tparam  maxOrder       max polynomial order
   * \tparam  Storage        select caching/non caching shape function set
   */
  template< class FunctionSpace, int maxOrder, template< class > class Storage >
  class ShapeFunctionSetStorage
  {
    typedef ShapeFunctionSetStorage< FunctionSpace, maxOrder, Storage > ThisType;

  public:
    typedef FunctionSpace FunctionSpaceType;

  private:
    typedef std::size_t size_type;

    static const int dimension = FunctionSpace::dimDomain;

    typedef MultiIndexSet< dimension, maxOrder > MultiIndexSetType;

    typedef ShapeFunctionSet< FunctionSpaceType, maxOrder > AnisotropicShapeFunctionSetType;

  public:
    typedef Dune::Fem::SelectCachingShapeFunctionSet< AnisotropicShapeFunctionSetType, Storage > ShapeFunctionSetType;
    typedef typename AnisotropicShapeFunctionSetType::MultiIndexType MultiIndexType;

  private:
    struct ShapeFunctionSetFactory
    {
      static Dune::GeometryType type ()
      {
        return Dune::GeometryType( typename Dune::GenericGeometry::CubeTopology< dimension >::type() );
      }

      static ShapeFunctionSetType *createObject ( const MultiIndexType &multiIndex )
      {
        return new ShapeFunctionSetType( type(), AnisotropicShapeFunctionSetType( multiIndex ) );
      }

      static void deleteObject ( ShapeFunctionSetType *shapeFunctionSet )
      {
        delete shapeFunctionSet;
      }
    };
    typedef ShapeFunctionSetFactory ShapeFunctionSetFactoryType;

    typedef Dune::Fem::SingletonList< MultiIndexType, ShapeFunctionSetType, ShapeFunctionSetFactoryType > SingletonProviderType;
   
  public:
    ShapeFunctionSetStorage ()
    : shapeFunctionSets_( MultiIndexSetType::size(), nullptr )
    {}

    ShapeFunctionSetStorage ( const ThisType &other )
    {
      typedef typename MultiIndexSetType::Iterator IteratorType;
      const IteratorType end = MultiIndexSetType::end();
      for( IteratorType it = MultiIndexSetType::begin(); it != end; ++it )
      {
        MultiIndexType multiIndex = *it;
        const size_type position = ThisType::position( multiIndex );

        ShapeFunctionSetType *shapeFunctionSet = other.shapeFunctionSets_[ position ];
        if( shapeFunctionSet )
        {
          shapeFunctionSets_[ position ] = &( SingletonProviderType::getObject( multiIndex ) );
        }
      }
    }

    ~ShapeFunctionSetStorage ()
    {
      typedef typename std::vector< ShapeFunctionSetType * >::iterator IteratorType;
      for( IteratorType it = shapeFunctionSets_.begin(); it != shapeFunctionSets_.end(); ++it )
      {
        ShapeFunctionSetType *shapeFunctionSet = *it;
        if( shapeFunctionSet )
        {
          SingletonProviderType::removeObject( *shapeFunctionSet );
          shapeFunctionSet = nullptr;
        }
      }
    }

    // add shape function set for given multi index to space
    bool insert ( const MultiIndexType &multiIndex )
    {
      const size_type position = ThisType::position( multiIndex );

      if( !shapeFunctionSets_[ position ] )
      {
        shapeFunctionSets_[ position ]
          = &( SingletonProviderType::getObject( multiIndex ) );
        assert( shapeFunctionSets_[ position ] );
        return true;
      }
      return false;
    }

    // return true, if shape function set for multi index has been inserted
    bool exists ( const MultiIndexType &multiIndex ) const
    {
      const size_type position = ThisType::position( multiIndex );
      return bool( shapeFunctionSets_[ position ] );
    }

    // return shape function set for multi index (do not test for null pointer)
    const ShapeFunctionSetType &operator[] ( const MultiIndexType &multiIndex ) const
    {
      const size_type position = ThisType::position( multiIndex );
      assert( shapeFunctionSets_[ position ] );
      return *shapeFunctionSets_[ position ];
    }

  private:
    // injective mapping from multi index to integer value
    size_type position ( const MultiIndexType &multiIndex ) const
    {
      std::size_t position = 0, factor = 1;
      for( std::size_t i = 0; i < dimension-1; ++i )
      {
        position += multiIndex[ i ]*factor;
        factor *= maxOrder+1;
      }
      position += multiIndex[ dimension-1 ]*factor;
      assert( position < MultiIndexSetType::size() );
      return position;
    }

    std::vector< ShapeFunctionSetType * > shapeFunctionSets_;
  };

} // namespace AnisotropicDG

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SHAPEFUNCTIONSETSTORAGE_HH
