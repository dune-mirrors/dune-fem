#ifndef DGADAPTIVELEAFGRIDPART_HH
#define DGADAPTIVELEAFGRIDPART_HH

//- Dune includes 
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/idbasedleafindexset.hh>
#include <dune/fem/storage/singletonlist.hh>

//- local includes 
#include "dgadaptiveleafindexset.hh"

namespace Dune
{

  // forward deklaration 
  template< class Grid >
  class DGAdaptiveLeafGridPart;

  //template< class Grid >
  //class DGAdaptiveLeafIndexSet;

  //! Type definitions for the DGAdaptiveLeafGridPart class
  template< class Grid >
  struct DGAdaptiveLeafGridPartTraits
  {
    //! default is to use DGAdaptiveLeafIndexSet
    template< class GridT, bool unStructured, bool isGood >
    struct GoodGridChooser
    {
      // choose the adative index based on hierarhic index set
      typedef DGAdaptiveLeafIndexSet< GridT > IndexSetType;
    };

    // the same for shitty unstructured grids 
    template< class GridT >
    struct GoodGridChooser< GridT, true, false >
    {
      // choose the leaf index set based on local ids 
      typedef IdBasedLeafIndexSet< GridT > IndexSetType;
    };

    // for structured grids with no hierarchic index set choose 
    // the DGAdaptiveLeafIndexSet based on the LeafIndexSet 
    template< class GridT >
    struct GoodGridChooser< GridT, false, false >
    {
      // choose the adative index based on leaf index set
      typedef DGAdaptiveLeafIndexSet< GridT > IndexSetType;
    };

    typedef Grid GridType;
    typedef DGAdaptiveLeafGridPart< GridType > GridPartType;

    // choose index set dependend on grid type  
    typedef GoodGridChooser
      < GridType, Capabilities :: IsUnstructured< GridType > :: v,
        Conversion< GridType, HasHierarchicIndexSet > :: exists >
      IndexSetChooserType;
                
    typedef typename IndexSetChooserType :: IndexSetType IndexSetType;

    static const PartitionIteratorType indexSetPartitionType = All_Partition;

    typedef typename GridType
      :: template Codim< 0 > :: Entity :: LeafIntersectionIterator
      IntersectionIteratorType;
    
    template< int codim >
    struct Codim
    {
      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef typename GridType
          :: template Codim< codim > :: template Partition< pitype > :: LeafIterator
          IteratorType;
      };
    };
        
    //! \brief is true if grid on this view only has conforming intersections 
    static const bool conforming = Capabilities :: isLeafwiseConforming< GridType > :: v;
  };


  /** @ingroup AdaptiveLeafGP
     \brief AdaptiveLeafGridPart with
     indexset only for codimension 0 entities.

      Special implementation of the AdaptiveLeafGridPart with
      an underlying index set is only defined for 
      entities with codimension 0 for use with
      the Dune::DiscontinuousGalerkinSpace.

      The underlying \ref DGAdaptiveLeafIndexSet "index set" is
      a singleton for each different grid. */
  template< class Grid >
  class DGAdaptiveLeafGridPart
  : public GridPartDefault< DGAdaptiveLeafGridPartTraits< Grid > >
  {
    typedef DGAdaptiveLeafGridPart< Grid > ThisType;
    typedef GridPartDefault< DGAdaptiveLeafGridPartTraits< Grid > >
      BaseType;

  public:
    //! Type definitions
    typedef DGAdaptiveLeafGridPartTraits< Grid > Traits;
    //! Grid implementation type
    typedef typename Traits :: GridType GridType;
    //! The leaf index set of the grid implementation
    typedef typename Traits :: IndexSetType IndexSetType;
    
    //! The corresponding IntersectionIterator 
    typedef typename Traits :: IntersectionIteratorType
      IntersectionIteratorType ;
    
    //! Struct providing types of the leaf iterators on codimension cd
    template< int codim >
    struct Codim
    : public BaseType :: template Codim< codim >
    {};

    typedef typename GridType::template Codim<0>::Entity EntityCodim0Type;

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

  public:
    //! Constructor
    explicit DGAdaptiveLeafGridPart( GridType &grid )
    : BaseType( grid, IndexSetProviderType :: getObject( &grid ) )
    {}
    //! Copy Constructor
    DGAdaptiveLeafGridPart( const ThisType &other )
    : BaseType( other.grid_, IndexSetProviderType :: getObject( &(other.grid()) ) )
    {}

    /** \brief Desctrutor removing index set. When only one reference is
        left, index set object is deleted. */
    inline ~DGAdaptiveLeafGridPart ()
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
    void communicate ( CommDataHandleIF< DataHandle, Data > &data,
                       InterfaceType iftype,
                       CommunicationDirection dir ) const
    {
      this->grid().communicate( data, iftype, dir );
    }
  };

}// end namespace Dune

#endif
