#ifndef DUNE_ADAPTIVELEAFGRIDPART_HH
#define DUNE_ADAPTIVELEAFGRIDPART_HH

#include <dune/common/typetraits.hh>

#include <dune/fem/misc/capabilities.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/gridpart/adaptiveleafindexset.hh>

#if DUNE_FEM_COMPATIBILITY
#include <dune/fem/space/dgspace/dgadaptiveleafgridpart.hh>
#endif

namespace Dune
{

  /** @addtogroup AdaptiveLeafGP
      GridPart for Dune::AdaptiveLeafIndexSet. 
      The underlying index set is
      singleton for each grid object.   
      Uses very efficient index sets specially
      designed for problems with constantly changing underlying grid.
  */
  /////////////////////////////////////////////////////////////////////////
  //
  //  --AdaptiveLeafIndexGridPart 
  //
  /////////////////////////////////////////////////////////////////////////

  // forward deklaration of grid part 
  template< class Grid, PartitionIteratorType idxpitype = All_Partition , bool onlyCodimensionZero = false >
  class AdaptiveLeafGridPart;


  //! Type definitions for the LeafGridPart class
  template< class Grid, PartitionIteratorType idxpitype , bool onlyCodimensionZero = false >
  class AdaptiveLeafGridPartTraits
  {
  public:
    //! type of the grid 
    typedef Grid GridType;

    //! type of the grid part , i.e. this type 
    typedef AdaptiveLeafGridPart< GridType, idxpitype, onlyCodimensionZero > GridPartType;

  protected:  
    // choose the AdaptiveIndexSet (based on the HierarchicIndexSet)
    // to be revised
    template < int dummy, bool onlyCodimZero > 
    struct AdaptiveLeafIndexSetChooser
    {
#ifdef USE_PARTITIONTYPED_INDEXSET
      static const PartitionIteratorType indexSetPartitionType = idxpitype;
#else
      static const PartitionIteratorType indexSetPartitionType = All_Partition;
#endif
      typedef AdaptiveLeafIndexSet< GridPartType > IndexSetType;
    };

#if ! DUNE_FEM_COMPATIBILITY // remove this #if in version 1.2
    template <int dummy> 
    struct AdaptiveLeafIndexSetChooser<dummy, true >
    {
#ifdef USE_PARTITIONTYPED_INDEXSET
      static const PartitionIteratorType indexSetPartitionType = idxpitype;
#else
      static const PartitionIteratorType indexSetPartitionType = All_Partition;
#endif
      typedef DGAdaptiveLeafIndexSet< GridPartType > IndexSetType;
    };
#endif

    // choose the LeafIndexSet
    struct LeafIndexSetChooser
    {
      static const PartitionIteratorType indexSetPartitionType = All_Partition;
      typedef WrappedLeafIndexSet< Grid > IndexSetType;
    };

    static const bool hasHierarchicIndexSet = Capabilities::hasHierarchicIndexSet< Grid >::v;
    typedef typename SelectType< hasHierarchicIndexSet, 
            AdaptiveLeafIndexSetChooser<-1, onlyCodimensionZero >, LeafIndexSetChooser>::Type
      IndexSetChooserType;

  public:  
    //! type of the index set 
    typedef typename IndexSetChooserType::IndexSetType IndexSetType;

    static const PartitionIteratorType indexSetPartitionType = IndexSetChooserType::indexSetPartitionType;

    typedef typename GridType::template Codim< 0 >::Entity::LeafIntersectionIterator
      IntersectionIteratorType;
    
    template< int cd >
    struct Codim
    {
      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef typename GridType
          :: template Codim< cd > :: template Partition< pitype > :: LeafIterator
          IteratorType;
      };
    };

    //! \brief is true if grid on this view only has conforming intersections
    static const bool conforming = Capabilities :: isLeafwiseConforming< GridType > :: v;
  };


  /** @ingroup AdaptiveLeafGP
      \brief A grid part with an index set specially
      designed for adaptive calculations.

      The underlying \ref AdaptiveLeafIndexSet "index set" is defined for
      entities of all codimensions. 
  */
  template< class Grid, PartitionIteratorType idxpitype , bool onlyCodimensionZero >
  class AdaptiveLeafGridPart
  : public GridPartDefault< AdaptiveLeafGridPartTraits< Grid, idxpitype, onlyCodimensionZero > >
  {
  public:  
    //! Type definitions
    typedef AdaptiveLeafGridPartTraits< Grid, idxpitype, onlyCodimensionZero > Traits;

  protected:  
    // type of this pointer 
    typedef AdaptiveLeafGridPart< Grid, idxpitype, onlyCodimensionZero > ThisType;

    // type of base class 
    typedef GridPartDefault< Traits > BaseType;

  public:
    //! Grid implementation type
    typedef typename Traits :: GridPartType GridPartType;
    //! Grid implementation type
    typedef typename Traits :: GridType GridType;
    //! The leaf index set of the grid implementation
    typedef typename Traits :: IndexSetType IndexSetType;

    //! The corresponding IntersectionIterator
    typedef typename Traits :: IntersectionIteratorType
      IntersectionIteratorType;

    typedef typename IntersectionIteratorType::Intersection IntersectionType;

    //! Struct providing types of the leaf iterators on codimension codim
    template< int codim >
    struct Codim
    : public BaseType :: template Codim< codim >
    {};

  private:
    struct IndexSetFactory
    {
      struct Key 
      {
        const GridPartType& gridPart_;
        const GridType& grid_;
        Key(const GridPartType& gridPart, const GridType& grid) 
         : gridPart_( gridPart ), grid_( grid ) 
        {}

        Key( const Key& other ) 
          : gridPart_( other.gridPart_ ) 
          , grid_( other.grid_ )
        {}
        bool operator ==( const Key& other ) const 
        {
          // compare grid pointers 
          return (&grid_) == (& other.grid_ );
        }
        const GridPartType& gridPart() const { return gridPart_; }
        const GridType& grid() const { return grid_; }
      };

      typedef IndexSetType ObjectType;
      typedef Key KeyType;

      inline static ObjectType *createObject ( const KeyType &key )
      {
        return new ObjectType( key.gridPart() );
      }

      inline static void deleteObject ( ObjectType *object )
      {
        delete object;
      }
    };

    typedef typename IndexSetFactory :: KeyType KeyType;
    typedef SingletonList
      < KeyType, IndexSetType, IndexSetFactory > IndexSetProviderType;

    // type of entity with codimension zero 
    typedef typename GridType :: template Codim< 0 > :: Entity EntityCodim0Type;

    // reference to index set 
    const IndexSetType& indexSet_;
  public:
    //! Constructor
    inline explicit AdaptiveLeafGridPart ( GridType &grid )
    : BaseType( grid )
    , indexSet_( IndexSetProviderType :: getObject( KeyType( *this, grid ) ) )
    {
    }

    //! Copy Constructor
    AdaptiveLeafGridPart ( const ThisType &other )
    : BaseType( other.grid_ )
    , indexSet_( IndexSetProviderType :: getObject( KeyType( *this, other.grid() ) ) )
    {
    }

    /** \brief Destrcutor removeing index set, if only one reference left, index set
        removed.  */
    inline ~AdaptiveLeafGridPart ()
    { 
      IndexSetProviderType :: removeObject( this->indexSet() );
    }

    //! Returns reference to index set of the underlying grid
    const IndexSetType &indexSet () const
    {
      return indexSet_;
    }

    //! Begin iterator on the leaf level
    template< int codim >
    typename Codim< codim > :: IteratorType
    begin () const
    {
      return BaseType :: template begin< codim >();
    }

    //! Begin iterator on the leaf level
    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition< pitype > :: IteratorType
    begin () const
    {
      return (*this).grid().template leafbegin< codim, pitype >();
    }

    //! Begin iterator on the leaf level
    template< int codim >
    typename Codim< codim > :: IteratorType
    end () const
    {
      return BaseType :: template end< codim >();
    }

    //! End iterator on the leaf level
    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim > :: template Partition< pitype > :: IteratorType
    end () const
    {
      return (*this).grid().template leafend< codim, pitype >();
    }

    //! ibegin of corresponding intersection iterator for given entity
    inline IntersectionIteratorType
    ibegin ( const EntityCodim0Type &entity ) const
    {
      return entity.ileafbegin();
    }
    
    //! iend of corresponding intersection iterator for given entity
    inline IntersectionIteratorType
    iend ( const EntityCodim0Type &entity ) const
    {
      return entity.ileafend();
    }

    int boundaryId ( const IntersectionType &intersection ) const
    {
      return intersection.boundaryId();
    }

    //! Returns maxlevel of the grid
    inline int level () const
    {
      return this->grid().maxLevel();
    }

    //! corresponding communication method for this grid part
    template< class DataHandle, class Data >
    inline void communicate ( CommDataHandleIF< DataHandle, Data > &data,
                              InterfaceType iftype,
                              CommunicationDirection dir ) const
    {
      this->grid().communicate( data, iftype, dir );
    }
  };

#if ! DUNE_FEM_COMPATIBILITY // remove this #if in version 1.2
  /** @ingroup AdaptiveLeafGP
      \brief A grid part with an index set specially
      designed for adaptive calculations.

      The underlying \ref DGAdaptiveLeafIndexSet "index set" is defined 
      only for codimension 0. 
  */
  template< class Grid, PartitionIteratorType idxpitype = All_Partition >
  class DGAdaptiveLeafGridPart
  : public AdaptiveLeafGridPart< Grid, idxpitype, true > 
  {
    typedef AdaptiveLeafGridPart< Grid, idxpitype, true > BaseType;
  public:  
    typedef typename BaseType :: GridType GridType;
    //! Constructor
    inline explicit DGAdaptiveLeafGridPart ( GridType &grid )
    : BaseType( grid )
    {
    }

    //! copy constructor 
    inline DGAdaptiveLeafGridPart ( const DGAdaptiveLeafGridPart& other )
    : BaseType( other )
    {
    }
  };
#endif

} // end namespace Dune

#endif
