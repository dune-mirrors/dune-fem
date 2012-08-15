#ifndef DUNE_FEM_ADAPTIVELEAFGRIDPART_HH
#define DUNE_FEM_ADAPTIVELEAFGRIDPART_HH

//- dune-common includes
#include <dune/common/typetraits.hh>

//- dune-fem includes
#include <dune/fem/gridpart/adaptiveleafindexset.hh>
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/defaultindexsets.hh>
#include <dune/fem/misc/capabilities.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/storage/singletonlist.hh>

namespace Dune
{

  namespace Fem
  {

    /////////////////////////////////////////////////////////////////////////
    //
    //  --AdaptiveLeafIndexGridPart 
    //
    /////////////////////////////////////////////////////////////////////////

    /** @ingroup AdaptiveLeafGP
        \brief A grid part with an index set specially
        designed for adaptive calculations.

        The underlying \ref AdaptiveLeafIndexSet "index set" is defined for
        entities of all codimensions. 
    */
    template< class TraitsImp >
    class AdaptiveGridPartBase
    : public GridPartDefault< TraitsImp >
    {
      typedef AdaptiveGridPartBase< TraitsImp > ThisType;
      typedef GridPartDefault< TraitsImp > BaseType;

    public:  
      //! Type definitions
      typedef TraitsImp Traits;

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
      typedef typename Codim< 0 > :: EntityType ElementType;

      typedef typename GridType :: template Partition < All_Partition > :: LeafGridView LeafGridView;
      // the leaf grid view 
      LeafGridView leafView_ ; 

      // reference to index set 
      const IndexSetType& indexSet_;

      // method to get DofManager instance to make sure the DofManager is delete after 
      // the index set provider 
      GridType &initDofManager(GridType &grid ) const 
      {
        DofManager< GridType > :: instance( grid );
        return grid ;
      }
    public:
      //! constructor
      explicit AdaptiveGridPartBase ( GridType &grid )
      : BaseType( initDofManager( grid ) ), // dofManager needs to be initialized before index set provider 
        leafView_( grid.leafView() ),
        indexSet_( IndexSetProviderType::getObject( KeyType( asImp(), grid ) ) )
      {}

      //! Copy Constructor
      AdaptiveGridPartBase ( const ThisType &other )
      : BaseType( other ),
        leafView_( other.leafView_ ),
        indexSet_( IndexSetProviderType::getObject( KeyType( asImp(), other.grid() ) ) )
      {}

      /** \brief Destrcutor removeing index set, if only one reference left, index set
          removed.  */
      ~AdaptiveGridPartBase ()
      { 
        IndexSetProviderType::removeObject( indexSet() );
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
        return begin< codim, InteriorBorder_Partition >();
      }

      //! Begin iterator on the leaf level
      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim > :: template Partition< pitype > :: IteratorType
      begin () const
      {
        return leafView_.template begin< codim, pitype >();
      }

      //! Begin iterator on the leaf level
      template< int codim >
      typename Codim< codim > :: IteratorType
      end () const
      {
        return end< codim, InteriorBorder_Partition >();
      }

      //! End iterator on the leaf level
      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim > :: template Partition< pitype > :: IteratorType
      end () const
      {
        return leafView_.template end< codim, pitype >();
      }

      //! ibegin of corresponding intersection iterator for given entity
      inline IntersectionIteratorType
      ibegin ( const ElementType &entity ) const
      {
        return entity.ileafbegin();
      }
      
      //! iend of corresponding intersection iterator for given entity
      inline IntersectionIteratorType
      iend ( const ElementType &entity ) const
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
    protected:
      const GridPartType& asImp() const 
      {
        return static_cast<const GridPartType &> (*this);
      }
    };

    /** @addtogroup AdaptiveLeafGP
        \brief GridPart for Dune::AdaptiveLeafIndexSet. 
        The underlying index set is
        singleton for each grid object.   
        Uses very efficient index sets specially
        designed for problems with constantly changing underlying grid.
    */
    template< class Grid, PartitionIteratorType idxpitype = All_Partition, bool onlyCodimensionZero = false >
    class AdaptiveLeafGridPart;

    //! Type definitions for the LeafGridPart class
    template< class Grid, PartitionIteratorType idxpitype , bool onlyCodimensionZero >
    class AdaptiveLeafGridPartTraits
    {
    public:
      //! type of the grid 
      typedef Grid GridType;

      //! type of the grid part , i.e. this type 
      typedef AdaptiveLeafGridPart< GridType, idxpitype, onlyCodimensionZero > GridPartType;

      /** \brief The type of the corresponding TwistUtility */
      typedef TwistUtility< GridType >  TwistUtilityType ;

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
        static const InterfaceType indexSetInterfaceType = All_All_Interface;
#endif
        typedef AdaptiveLeafIndexSet< GridPartType > IndexSetType;
      };

      template <int dummy> 
      struct AdaptiveLeafIndexSetChooser<dummy, true >
      {
#ifdef USE_PARTITIONTYPED_INDEXSET
        static const PartitionIteratorType indexSetPartitionType = idxpitype;
#else
        static const PartitionIteratorType indexSetPartitionType = All_Partition;
        static const InterfaceType indexSetInterfaceType = All_All_Interface;
#endif
        typedef DGAdaptiveLeafIndexSet< GridPartType > IndexSetType;
      };

      // also for Cartesian grids (e.g. YaspGrid) use adaptive leaf index set in parallel 
      typedef AdaptiveLeafIndexSetChooser<-1, onlyCodimensionZero > IndexSetChooserType;

    public:  
      //! type of the index set 
      typedef typename IndexSetChooserType::IndexSetType IndexSetType;

      static const PartitionIteratorType indexSetPartitionType = IndexSetChooserType::indexSetPartitionType;
      static const InterfaceType indexSetInterfaceType = IndexSetChooserType::indexSetInterfaceType;

      typedef typename GridType::template Codim< 0 >::Entity::LeafIntersectionIterator
        IntersectionIteratorType;
      
      template< int codim >
      struct Codim
      {
        typedef typename GridType::template Codim< codim >::Geometry      GeometryType;
        typedef typename GridType::template Codim< codim >::LocalGeometry LocalGeometryType;

        typedef typename GridType::template Codim< codim >::EntityPointer EntityPointerType;
        typedef typename GridType::template Codim< codim >::Entity        EntityType;

        typedef typename GridType::template Codim< codim >::EntitySeed    EntitySeedType;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef typename GridType::template Codim< codim >::template Partition< pitype >::LeafIterator IteratorType;
        };
      };

      //! \brief is true if grid on this view only has conforming intersections
      static const bool conforming = Capabilities :: isLeafwiseConforming< GridType > :: v;
    };

    template< class Grid, PartitionIteratorType idxpitype , bool onlyCodimensionZero >
    class AdaptiveLeafGridPart
      : public AdaptiveGridPartBase< AdaptiveLeafGridPartTraits< Grid, idxpitype, onlyCodimensionZero > >
    {
      typedef AdaptiveGridPartBase< AdaptiveLeafGridPartTraits< Grid, idxpitype, onlyCodimensionZero > > BaseType;
    public:  
      typedef typename BaseType :: GridType GridType;
      //! Constructor
      inline explicit AdaptiveLeafGridPart ( GridType &grid )
      : BaseType( grid )
      {
      }

      //! copy constructor 
      inline AdaptiveLeafGridPart ( const AdaptiveLeafGridPart& other )
      : BaseType( other )
      {
      }
    };

    /** @ingroup AdaptiveLeafGP
        \brief A grid part with an index set specially
        designed for adaptive calculations.

        The underlying \ref DGAdaptiveLeafIndexSet "index set" is defined 
        only for codimension 0. 
    */
    template< class Grid, PartitionIteratorType idxpitype = All_Partition >
    class DGAdaptiveLeafGridPart
    : public AdaptiveGridPartBase< AdaptiveLeafGridPartTraits< Grid, idxpitype, true > > 
    {
      typedef AdaptiveGridPartBase< AdaptiveLeafGridPartTraits< Grid, idxpitype, true > > BaseType;
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

    template< class Grid, PartitionIteratorType idxpitype = All_Partition >
    class IntersectionAdaptiveLeafGridPart ;

    //! Type definitions for the LeafGridPart class
    template< class Grid, PartitionIteratorType idxpitype >
    class IntersectionAdaptiveLeafGridPartTraits : public AdaptiveLeafGridPartTraits< Grid, idxpitype, false>
    {
    public:
      //! type of the grid part , i.e. this type 
      typedef IntersectionAdaptiveLeafGridPart< Grid, idxpitype > GridPartType;

      //! type of the index set 
      typedef IntersectionAdaptiveLeafIndexSet< GridPartType > IndexSetType;
    };

    /** @ingroup AdaptiveLeafGP
        \brief A grid part with an index set specially
        designed for adaptive calculations.

        The underlying \ref DGAdaptiveLeafIndexSet "index set" is defined 
        only for codimension 0. 
    */
    template< class Grid, PartitionIteratorType idxpitype >
    class IntersectionAdaptiveLeafGridPart
    : public AdaptiveGridPartBase< IntersectionAdaptiveLeafGridPartTraits< Grid, idxpitype > > 
    {
      typedef AdaptiveGridPartBase< IntersectionAdaptiveLeafGridPartTraits< Grid, idxpitype > > BaseType;
    public:  
      typedef typename BaseType :: GridType GridType;
      //! Constructor
      inline explicit IntersectionAdaptiveLeafGridPart( GridType &grid )
      : BaseType( grid )
      {
      }

      //! copy constructor 
      inline IntersectionAdaptiveLeafGridPart( const IntersectionAdaptiveLeafGridPart& other )
      : BaseType( other )
      {
      }
    };



    // Capabilities
    // ------------

    namespace GridPartCapabilities
    {

      // Capabilities for AdaptiveLeafGridPart
      // -------------------------------------

      template< class Grid, PartitionIteratorType idxpitype, bool onlyCodimensionZero >
      struct hasGrid< AdaptiveLeafGridPart< Grid, idxpitype, onlyCodimensionZero > >
      {
        static const bool v = true;
      };

      template< class Grid, PartitionIteratorType idxpitype, bool onlyCodimensionZero >
      struct hasSingleGeometryType< AdaptiveLeafGridPart< Grid, idxpitype, onlyCodimensionZero > >
      {
        static const bool v = Dune::Capabilities::hasSingleGeometryType< Grid >::v;
        static const unsigned int topologyId
          = Dune::Capabilities::hasSingleGeometryType< Grid >::topologyId;;
      };

      template< class Grid, PartitionIteratorType idxpitype, bool onlyCodimensionZero >
      struct isCartesian< AdaptiveLeafGridPart< Grid, idxpitype, onlyCodimensionZero > >
      {
        static const bool v = Dune::Capabilities::isCartesian< Grid >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype, bool onlyCodimensionZero, int codim  >
      struct hasEntity< AdaptiveLeafGridPart< Grid, idxpitype, onlyCodimensionZero >, codim >
      {
        static const bool v = Dune::Capabilities::hasEntity< Grid, codim >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype, bool onlyCodimensionZero >
      struct isParallel< AdaptiveLeafGridPart< Grid, idxpitype, onlyCodimensionZero > >
      {
        static const bool v = Dune::Capabilities::isParallel< Grid >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype, bool onlyCodimensionZero, int codim >
      struct canCommunicate< AdaptiveLeafGridPart< Grid, idxpitype, onlyCodimensionZero >, codim >
      {
        static const bool v = Dune::Capabilities::canCommunicate< Grid, codim >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype, bool onlyCodimensionZero >
      struct isConforming< AdaptiveLeafGridPart< Grid, idxpitype, onlyCodimensionZero > >
      {
        static const bool v = Dune::Capabilities::isLeafwiseConforming< Grid >::v;
      };


      // Capabilities for DGAdaptiveLeafGridPart
      // ---------------------------------------

      template< class Grid, PartitionIteratorType idxpitype >
      struct hasGrid< DGAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = true;
      };

      template< class Grid, PartitionIteratorType idxpitype >
      struct hasSingleGeometryType< DGAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = Dune::Capabilities::hasSingleGeometryType< Grid >::v;
        static const unsigned int topologyId
          = Dune::Capabilities::hasSingleGeometryType< Grid >::topologyId;;
      };

      template< class Grid, PartitionIteratorType idxpitype >
      struct isCartesian< DGAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = Dune::Capabilities::isCartesian< Grid >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype, int codim  >
      struct hasEntity< DGAdaptiveLeafGridPart< Grid, idxpitype >, codim >
      {
        static const bool v = Dune::Capabilities::hasEntity< Grid, codim >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype >
      struct isParallel< DGAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = Dune::Capabilities::isParallel< Grid >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype, int codim >
      struct canCommunicate< DGAdaptiveLeafGridPart< Grid, idxpitype >, codim >
      {
        static const bool v = Dune::Capabilities::canCommunicate< Grid, codim >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype >
      struct isConforming< DGAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = Dune::Capabilities::isLeafwiseConforming< Grid >::v;
      };


      // Capbilities for IntersectionAdaptiveLeafGridPart
      // ------------------------------------------------

      template< class Grid, PartitionIteratorType idxpitype >
      struct hasGrid< IntersectionAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = true;
      };

      template< class Grid, PartitionIteratorType idxpitype >
      struct hasSingleGeometryType< IntersectionAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = Dune::Capabilities::hasSingleGeometryType< Grid >::v;
        static const unsigned int topologyId
          = Dune::Capabilities::hasSingleGeometryType< Grid >::topologyId;;
      };

      template< class Grid, PartitionIteratorType idxpitype >
      struct isCartesian< IntersectionAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = Dune::Capabilities::isCartesian< Grid >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype, int codim  >
      struct hasEntity< IntersectionAdaptiveLeafGridPart< Grid, idxpitype >, codim >
      {
        static const bool v = Dune::Capabilities::hasEntity< Grid, codim >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype >
      struct isParallel< IntersectionAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = Dune::Capabilities::isParallel< Grid >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype, int codim >
      struct canCommunicate< IntersectionAdaptiveLeafGridPart< Grid, idxpitype >, codim >
      {
        static const bool v = Dune::Capabilities::canCommunicate< Grid, codim >::v;
      };

      template< class Grid, PartitionIteratorType idxpitype >
      struct isConforming< IntersectionAdaptiveLeafGridPart< Grid, idxpitype > >
      {
        static const bool v = Dune::Capabilities::isLeafwiseConforming< Grid >::v;
      };

    } // namespace GridPartCapabilities

  } // namespace Fem

#if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: AdaptiveLeafGridPart ;
  using Fem :: DGAdaptiveLeafGridPart ;
  using Fem :: IntersectionAdaptiveLeafGridPart ;
#endif // DUNE_FEM_COMPATIBILITY

} // namespace Dune

#endif // #ifndef DUNE_FEM_ADAPTIVELEAFGRIDPART_HH
