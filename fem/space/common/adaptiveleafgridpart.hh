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
  template< class Grid, PartitionIteratorType pitype >
  class AdaptiveLeafGridPart;

  template< class Grid >
  class AdaptiveLeafIndexSet;

  //! Type definitions for the LeafGridPart class
  template< class Grid, PartitionIteratorType pitype >
  struct AdaptiveLeafGridPartTraits
  {
    //! type of the grid 
    typedef Grid GridType;

    //! default is to use DGAdaptiveLeafIndexSet
    template< class GridT, bool isGood >
    struct GoodGridChooser
    {
      // choose the adative index based on hierarhic index set
      typedef AdaptiveLeafIndexSet< GridT > IndexSetType;
    };

    // the same for shitty grids 
    template< class GridT >
    struct GoodGridChooser< GridT, false >
    {
      // the grids leaf index set wrapper for good 
      typedef WrappedLeafIndexSet< GridT > IndexSetType;
    };

    //! type of the grid part , i.e. this type 
    typedef AdaptiveLeafGridPart< GridType, pitype > GridPartType;

    // choose index set dependend on grid type  
    typedef GoodGridChooser
      < GridType, Conversion< GridType, HasHierarchicIndexSet > :: exists >
      IndexSetChooserType;
                
    //! type of the index set 
    typedef typename IndexSetChooserType :: IndexSetType IndexSetType;

    typedef typename GridType
      :: template Codim< 0 > :: Entity :: LeafIntersectionIterator
      IntersectionIteratorType;
    
    template< int cd >
    struct Codim
    {
      typedef typename GridType
        :: template Codim< cd > :: template Partition< pitype > :: LeafIterator
        IteratorType;
    };

    //! \brief is true if grid on this view only has conforming intersections 
    enum { conforming = Capabilities :: isLeafwiseConforming< GridType > :: v };
  };


  /** @ingroup AdaptiveLeafGP
      \brief A grid part with an index set specially
      designed for adaptive calculations.

      The underlying \ref AdaptiveLeafIndexSet "index set" is defined for
      entities of all codimensions. 
  */
  template< class Grid, PartitionIteratorType pitype = Interior_Partition >
  class AdaptiveLeafGridPart
  : public GridPartDefault< AdaptiveLeafGridPartTraits< Grid, pitype > >
  {
    typedef AdaptiveLeafGridPart< Grid, pitype > ThisType;
    typedef GridPartDefault< AdaptiveLeafGridPartTraits< Grid, pitype > >
      BaseType;

  public:
    //! Type definitions
    typedef AdaptiveLeafGridPartTraits< Grid, pitype > Traits;

    //! Grid implementation type
    typedef typename Traits :: GridType GridType;
    //! The leaf index set of the grid implementation
    typedef typename Traits :: IndexSetType IndexSetType;

    //! The corresponding IntersectionIterator
    typedef typename Traits :: IntersectionIteratorType
      IntersectionIteratorType;

    //! Struct providing types of the leaf iterators on codimension cd
    template< int cd >
    struct Codim
    {
      typedef typename Traits :: template Codim< cd > :: IteratorType
        IteratorType;
    };

    //! \brief is true if grid on this view only has conforming intersections 
    enum { conforming = Traits :: conforming };

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

    /** \brief Destrcutor removeing index set, if only one reference left, index set
        removed.  */
    inline ~AdaptiveLeafGridPart ()
    { 
      IndexSetProviderType :: removeObject( this->indexSet() );
    }

    //! Begin iterator on the leaf level
    template< int cd >
    inline typename Codim< cd > :: IteratorType begin () const
    {
      return this->grid().template leafbegin< cd, pitype >();
    }

    //! End iterator on the leaf level
    template< int cd >
    inline typename Codim< cd > :: IteratorType end () const
    {
      return this->grid().template leafend< cd, pitype >();
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
