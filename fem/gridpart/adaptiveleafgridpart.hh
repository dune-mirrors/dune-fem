#ifndef DUNE_ADAPTIVELEAFGRIDPART_HH
#define DUNE_ADAPTIVELEAFGRIDPART_HH

//- Dune includes 
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/storage/singletonlist.hh>

//- local includes 
#include "adaptiveleafindexset.hh"

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
  template< class Grid, PartitionIteratorType idxpitype = All_Partition >
  class AdaptiveLeafGridPart;

  //template< class Grid, PartitionIteratorType pitype >
  //class AdaptiveLeafIndexSet;



  //! Type definitions for the LeafGridPart class
  template< class Grid, PartitionIteratorType idxpitype >
  struct AdaptiveLeafGridPartTraits
  {
    //! type of the grid 
    typedef Grid GridType;

    //! default is to use DGAdaptiveLeafIndexSet
    template< class GridT, bool isGood >
    struct GoodGridChooser
    {
      // choose the adative index based on hierarhic index setq
      // to be revised 
#ifdef USE_PARTITIONTYPED_INDEXSET
      static const PartitionIteratorType indexSetPartitionType = idxpitype;
#else
      static const PartitionIteratorType indexSetPartitionType = All_Partition;
#endif
      typedef AdaptiveLeafIndexSet< GridT, indexSetPartitionType > IndexSetType;
    };

    // the same for shitty grids 
    template< class GridT >
    struct GoodGridChooser< GridT, false >
    {
      static const PartitionIteratorType indexSetPartitionType = All_Partition;
      // the grids leaf index set wrapper for good 
      typedef WrappedLeafIndexSet< GridT > IndexSetType;
    };

    //! type of the grid part , i.e. this type 
    typedef AdaptiveLeafGridPart< GridType, idxpitype > GridPartType;

    // choose index set dependend on grid type  
    typedef GoodGridChooser
      < GridType, Conversion< GridType, HasHierarchicIndexSet > :: exists >
      IndexSetChooserType;
                
    //! type of the index set 
    typedef typename IndexSetChooserType :: IndexSetType IndexSetType;

    static const PartitionIteratorType indexSetPartitionType
      = IndexSetChooserType :: indexSetPartitionType;

    typedef typename GridType
      :: template Codim< 0 > :: Entity :: LeafIntersectionIterator
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
  template< class Grid, PartitionIteratorType idxpitype >
  class AdaptiveLeafGridPart
  : public GridPartDefault< AdaptiveLeafGridPartTraits< Grid, idxpitype > >
  {
    typedef AdaptiveLeafGridPart< Grid, idxpitype > ThisType;
    typedef GridPartDefault< AdaptiveLeafGridPartTraits< Grid, idxpitype > >
      BaseType;

  public:
    //! Type definitions
    typedef AdaptiveLeafGridPartTraits< Grid, idxpitype > Traits;

    //! Grid implementation type
    typedef typename Traits :: GridType GridType;
    //! The leaf index set of the grid implementation
    typedef typename Traits :: IndexSetType IndexSetType;

    //! The corresponding IntersectionIterator
    typedef typename Traits :: IntersectionIteratorType
      IntersectionIteratorType;

    //! Struct providing types of the leaf iterators on codimension codim
    template< int codim >
    struct Codim
    : public BaseType :: template Codim< codim >
    {};

  private:
    struct IndexSetFactory
    {
      typedef IndexSetType ObjectType;
      typedef const GridType *KeyType;

      inline static ObjectType *createObject ( const KeyType &key )
      {
        return new ObjectType( *key );
      }

      inline static void deleteObject ( ObjectType *object )
      {
        delete object;
      }
    };

    typedef SingletonList
      < typename IndexSetFactory :: KeyType, IndexSetType, IndexSetFactory >
      IndexSetProviderType;

    typedef typename GridType :: template Codim< 0 > :: Entity EntityCodim0Type;

  public:
    //! Constructor
    inline explicit AdaptiveLeafGridPart ( GridType &grid )
    : BaseType( grid, IndexSetProviderType :: getObject( &grid ) )
    {}
    //! Copy Constructor
    AdaptiveLeafGridPart ( const ThisType &other )
    : BaseType( other.grid_, IndexSetProviderType :: getObject( &(other.grid()) ) )
    {}

    /** \brief Destrcutor removeing index set, if only one reference left, index set
        removed.  */
    inline ~AdaptiveLeafGridPart ()
    { 
      IndexSetProviderType :: removeObject( this->indexSet() );
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

} // end namespace Dune

#endif
