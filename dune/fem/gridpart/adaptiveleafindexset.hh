#ifndef DUNE_FEM_ADAPTIVELEAFINDEXSET_HH
#define DUNE_FEM_ADAPTIVELEAFINDEXSET_HH

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <string>
#include <vector>
#include <memory>
#include <type_traits>

#include <dune/fem/common/forloop.hh>

#include <dune/fem/gridpart/codimindexset.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/gridpart/common/persistentindexset.hh>
#include <dune/fem/io/file/iointerface.hh>
#include <dune/fem/version.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal forward declations
    // ---------------------------

    template < class GridPartImp >
    class AdaptiveLeafIndexSet;
    template < class GridPartImp >
    class IntersectionAdaptiveLeafIndexSet;
    template < class GridPartImp >
    class DGAdaptiveLeafIndexSet;




    /////////////////////////////////////////////////////////////////////////
    //
    //  --AdaptiveIndexSetBaseTraits
    //
    /////////////////////////////////////////////////////////////////////////

    template< class GridPart, class IndexSet >
    class AdaptiveIndexSetBaseTraits
    {
    public:
      // index set type derived from AdaptiveIndexSetBase
      typedef IndexSet IndexSetType;

      // grid part type
      typedef GridPart GridPartType;
      // grid type
      typedef typename GridPartType :: GridType GridType;

      // grid (part) dimension
      static const int dimension = GridPartType :: dimension;

      template< int codim >
      struct Codim
      {
        // entity type
        typedef typename GridPartType :: template Codim< codim > :: EntityType Entity;
      };

      // type of codim index set
      typedef CodimIndexSet< GridType > CodimIndexSetType;
      // index type used
      typedef typename CodimIndexSetType :: IndexType IndexType;

      // container of geometry types
      typedef std::vector< GeometryType > Types;
    };



    /////////////////////////////////////////////////////////////////////////
    //
    //  --AdaptiveIndexSetBase
    //
    /////////////////////////////////////////////////////////////////////////

    /** \class AdaptiveIndexSetBase
     *  \brief consecutive, persistent index set for the leaf level based on the
     *         grid's hierarchy index set
     *
     *  This index set generates a consecutive leaf index out of the unique global
     *  index of each entity. It can be used instead of the default grid index sets
     *  and can be generated for each grid implementation.
     *
     *  \note The base implementation can support either only one codimension or all
     *  codimensions of the grid.
     */
    template <class TraitsImp >
    class AdaptiveIndexSetBase
    : public PersistentAdaptiveIndexSet< TraitsImp >
    {
      typedef AdaptiveIndexSetBase< TraitsImp > ThisType;
      typedef PersistentAdaptiveIndexSet< TraitsImp > BaseType;

    protected:
      typedef typename TraitsImp :: GridPartType GridPartType;
      typedef typename GridPartType :: GridType GridType;

      typedef typename TraitsImp :: CodimIndexSetType  CodimIndexSetType ;

      typedef typename GridType::template Codim< 0 >::Entity GridElementType;

    public:
      //! \copydoc Dune::Fem::IndexSet::dimension */
      static const int dimension = BaseType::dimension;

      //! number of supported codimensions
      static const int numCodimensions = TraitsImp :: numCodimensions ;

      //! intersection codimension (numCodim-1 if enabled, otherwise -1)
      static const int intersectionCodimension = TraitsImp :: intersectionCodimension ;

      //! true if only one geometry type is available
      static const bool hasSingleGeometryType = Dune::Capabilities::hasSingleGeometryType< GridType > :: v ;

      //! \copydoc Dune::Fem::IndexSet::IndexType */
      typedef typename BaseType :: IndexType IndexType;

      //! \copydoc Dune::Fem::IndexSet::Types */
      typedef typename BaseType :: Types Types;

      //! type of codimension 0 Entity
      typedef typename BaseType :: template Codim< 0 > :: Entity ElementType;

      //! type of intersection iterator
      typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;

      //! type of intersections
      typedef typename GridPartType :: IntersectionType   IntersectionType;


    private:

      template< int codim , bool gridHasCodim >
      struct CountElementsBase
      {
        static void apply ( const ThisType &indexSet, const GeometryType &type, IndexType& count )
        {
          if( type.dim() == dimension - codim )
            count = indexSet.template countElements< codim >( type, std::integral_constant<bool,true>() );
        }
      };

      template< int codim >
      struct CountElementsBase< codim, false >
      {
        static void apply ( const ThisType &indexSet, const GeometryType &type, IndexType& count )
        {
          if( type.dim() == dimension - codim )
            count = indexSet.template countElements< codim >( type, std::integral_constant<bool,false>() );
        }
      };

      template< int codim >
      struct CountElements
        : public CountElementsBase< codim, Dune::Fem::GridPartCapabilities::hasEntity< GridPartType, codim > :: v >
      {
      };


      template< int codim >
      struct InsertSubEntities
      {
        static void apply ( const ThisType &indexSet, const GridElementType &element )
        {
          // if codimension is not available return
          if( ! indexSet.codimAvailable( codim ) ) return ;

          // if codimension is not used return
          if( !indexSet.codimUsed( codim ) ) return;

          CodimIndexSetType &codimSet = indexSet.codimLeafSet( codim );

          const int count = element.subEntities( codim );
          for( int i = 0; i < count; ++i )
          {
            codimSet.insertSubEntity( element, i );
          }
        }
      };

      template< int codim , bool gridHasCodim >
      struct InsertGhostSubEntitiesBase
      {
        static void apply ( ThisType &indexSet, const GridElementType &element,
                            const bool skipGhosts )
        {
          // if codimension is not available return
          if( ! indexSet.codimAvailable( codim ) ) return ;

          // if codimension is not used return
          if( !indexSet.codimUsed( codim ) ) return;

          CodimIndexSetType &codimSet = indexSet.codimLeafSet( codim );

          for( unsigned int i = 0; i < element.subEntities( codim ); ++i )
          {
            if( !skipGhosts || (element.partitionType() != GhostEntity) )
              codimSet.insertGhost( element.template subEntity< codim >( i ) );
          }
        }
      };

      template< int codim >
      struct InsertGhostSubEntitiesBase< codim, false >
      {
        static void apply ( ThisType &indexSet, const GridElementType &element,
                            const bool skipGhosts )
        {}
      };

      template< int codim >
      struct InsertGhostSubEntities
        : public InsertGhostSubEntitiesBase< codim, Dune::Capabilities::hasEntity < GridType, codim > :: v >
      {
      };

      template< int codim , bool gridHasCodim >
      struct CallSetUpCodimSetBase
      {
        static void apply ( const int cd, const ThisType &indexSet )
        {
          // if codimension is not available return
          if( ! indexSet.codimAvailable( codim ) ) return ;

          if( cd == codim && ! indexSet.codimUsed( cd ) )
            indexSet.template setupCodimSet< codim >(std::integral_constant<bool,true>());
        }
      };

      template< int codim >
      struct CallSetUpCodimSetBase< codim, false >
      {
        static void apply ( const int cd, const ThisType &indexSet )
        {
          if( cd == codim && ! indexSet.codimUsed( cd ) )
            indexSet.template setupCodimSet< codim >(std::integral_constant<bool,false>());
        }
      };

      template< int codim >
      struct CallSetUpCodimSet
        : public CallSetUpCodimSetBase< codim, Dune::Capabilities::hasEntity < GridType, codim > :: v >
      {
      };


      /////////////////////////////////////////////////////
      //  subentity extractor
      /////////////////////////////////////////////////////

      template < int codim, bool gridHasCodim >
      struct GetSubEntityBase
      {
        typedef typename GridPartType :: template Codim< codim > :: EntityType  Entity;
        static Entity subEntity( const ElementType& element, const int subEn )
        {
          return element.template subEntity< codim > ( subEn );
        }
      };

      template < int codim >
      struct GetSubEntityBase< codim, false >
      {
        typedef typename GridPartType :: template Codim< 0 > :: EntityType  Entity;
        static Entity subEntity( const ElementType& element, const int subEn )
        {
          DUNE_THROW(NotImplemented,"stupid grid without entities of codim 1 used");
        }
      };

      struct GetFaceEntity
        : public GetSubEntityBase< 1, Dune::Capabilities::hasEntity < GridType, 1 > :: v  >
      {
      };

      //! type of codimension 1 entity
      typedef typename GetFaceEntity :: Entity FaceType;

      //! is true if grid is a Cartesian, non-adaptive grid (i.e. YaspGrid, SPGrid)
      enum { CartesianNonAdaptiveGrid =  Dune::Capabilities::isCartesian<GridType>::v &&
                                        ! Capabilities::isLocallyAdaptive<GridType>::v };

      // my type, to be revised
      enum { myType = ( numCodimensions == 1 ) ? ( (CartesianNonAdaptiveGrid) ? -1 : 665 ) : 6 };

      // max num of codimension (to avoid compiler warnings)
      enum { maxNumCodimension = ((dimension + 1) > numCodimensions) ? dimension + 2 : numCodimensions+1 };

      //! default partition iterator type
      static const PartitionIteratorType pitype = GridPartType :: indexSetPartitionType ;

      // grid part (for iterator access, no index set)
      const GridPartType& gridPart_;

      // pointer storage in case index set is created by passing a grid
      std::unique_ptr< const GridPartType > gridPartPtr_;

      // Codimension leaf index sets
      mutable std::unique_ptr< CodimIndexSetType > codimLeafSet_[ numCodimensions ];

      // vector holding geometry types
      std::vector< std::vector< GeometryType > > geomTypes_;

      // actual sequence number
      int sequence_;

      //! flag is tru if set is in compressed status
      mutable bool compressed_;

    protected:
      using BaseType::grid_;
      using BaseType::dofManager_;

      // return true if codim is supported
      bool codimAvailable( const int codim ) const
      {
        return (codim < numCodimensions && codim >= 0);
      }

      // return true if indices for this codim exist
      bool codimUsed( const int codim ) const
      {
        return codimAvailable( codim ) && codimLeafSet_[ codim ] ;
      }

      CodimIndexSetType& codimLeafSet( const int codim ) const
      {
        assert( codimUsed( codim ) );
        return *codimLeafSet_[ codim ];
      }

    public:
      void requestCodimensions ( const std::vector< int >& codimensions ) const
      {
        // enable requested codimensions and rebuild index set
        for( const auto& codim : codimensions )
        {
          // start loop from 1 since codim 0 is always created
          Fem::ForLoop< CallSetUpCodimSet, 1, dimension >::apply( codim, *this );
        }
      }

      //! Constructor
      AdaptiveIndexSetBase ( const GridType* grid )
        : AdaptiveIndexSetBase( *(new GridPartType( const_cast< GridType& > (*grid), typename GridPartType::NoIndexSetType() )) )
      {
        // store pointer to avoid memory leaks
        gridPartPtr_.reset( &gridPart_ );
      }

      //! Constructor
      AdaptiveIndexSetBase ( const GridPartType& gridPart )
        : BaseType( gridPart.grid() )
        , gridPart_( gridPart )
        , sequence_( dofManager_.sequence() )
        , compressed_(true) // at start the set is compressed
      {
        // codim 0 is used by default
        codimLeafSet_[ 0 ].reset( new CodimIndexSetType( grid_, 0 ) );

        // set the codim of each Codim Set.
        for(int codim = 0; codim < numCodimensions; ++codim )
        {
          if( codim == intersectionCodimension )
            codimLeafSet_[ codim ].reset( new CodimIndexSetType( grid_, 1 ) );
        }

        /// get geometry types (not working for hybrid grids, like to whole set itself)
        {
          // get level-0 view, this is already used in GridPtr (DFG parser)
          typedef typename GridType :: LevelGridView MacroViewType;
          MacroViewType macroView = grid_.levelGridView( 0 );
          const typename MacroViewType :: IndexSet& indexSet = macroView.indexSet();

          // resize vector of geometry types
          geomTypes_.resize( dimension+1 );
          for(int codim=0; codim <= dimension; ++codim )
          {
            const int size = indexSet.types( codim ).size();
            // copy geometry types
            geomTypes_[ codim ].resize( size );
            std::copy_n( indexSet.types( codim ).begin(), size, geomTypes_[ codim ].begin() );
          }
        }

        // build index set
        setupIndexSet();
      }

      //! return type of index set, for GrapeDataIO
      int type () const
      {
        return myType;
      }

      //! return name of index set
      virtual std::string name () const
      {
        return "AdaptiveIndexSetBase";
      }

      //****************************************************************
      //
      //  INTERFACE METHODS for DUNE INDEX SETS
      //  --size
      //
      //****************************************************************
      //! \copydoc Dune::Fem::IndexSet::size */
      IndexType size ( GeometryType type ) const
      {
        const int codim = dimension - type.dim();

        // true if only one geometry type is present
        const bool onlySingleGeometryType = hasSingleGeometryType || ( geomTypes( codim ).size() == 1 ) ;
        // use size of codim index set if possible
        if( codimAvailable( codim ) && onlySingleGeometryType )
        {
          if( codimUsed( codim ) )
            return type == geomTypes( codim )[ 0 ] ? codimLeafSet( codim ).size() : 0;
        }

        // count entities for given geometry type
        IndexType count = 0 ;
        Fem::ForLoop< CountElements, 0, dimension > :: apply( *this, type, count );
        return count;
      }

      //! \copydoc Dune::Fem::IndexSet::size */
      IndexType size ( int codim ) const
      {
        assert( codim < numCodimensions );
        if( intersectionCodimension > 0 && codim == intersectionCodimension )
        {
          return codimLeafSet( codim ).size();
        }

        // count size for all geometry types
        IndexType count = 0 ;
        const size_t types = geomTypes( codim ).size();
        for( size_t i=0; i<types; ++i )
        {
          count += size( geomTypes( codim )[ i ] );
        }
        return count ;
      }

      //! \copydoc Dune::Fem::IndexSet::geomTypes */
      const std::vector <GeometryType> & geomTypes (const int codim) const
      {
        assert( codim >= 0 && codim < int(geomTypes_.size()) );
        return geomTypes_[ codim ];
      }

      //! \copydoc Dune::Fem::IndexSet::types */
      Types types( const int codim ) const
      {
        return geomTypes( codim );
      }

      //! \copydoc Dune::Fem::IndexSet::contains */
      template <class EntityType>
      bool contains (const EntityType & en) const
      {
        enum { codim = EntityType::codimension };
        if( codimAvailable( codim ) )
        {
          assert( codimUsed( codim ) );
          return codimLeafSet( codim ).exists( gridEntity( en ) );
        }
        else
          return false;
      }

      //****************************************************************
      //
      //  METHODS for Adaptation with DofManger
      //
      //****************************************************************

      //! \copydoc Dune::Fem::ConsecutiveIndexSet::insertEntity */
      void insertEntity( const GridElementType &entity )
      {
        // here we have to add the support of higher codims
        resizeVectors();
        insertIndex( entity );
      }

      //! \copydoc Dune::Fem::ConsecutiveIndexSet::removeEntity */
      void removeEntity( const GridElementType &entity )
      {
        removeIndex( entity );
      }

      //! reallocate the vector for new size
      void resizeVectors ();

      //! \copydoc Dune::Fem::ConsecutiveIndexSet::resize */
      void resize ()
      {
        resizeVectors();

    #if HAVE_MPI
        if( CartesianNonAdaptiveGrid  &&
            grid_.comm().size() > 1 )
        {
          // only done for structured grids
          clear();

          // this should only be the case of YaspGrid
          markAllBelowOld<Interior_Partition>();
          if( pitype > Interior_Partition )
          {
            markAllBelowOld< pitype >();
          }
          compressed_ = true;
        }
        else
    #endif
        {
          // use a hierarchic walk to mark new elements
          markAllBelowOld< pitype > ();

    #if HAVE_MPI
          // only if ghost are really supported
          if( pitype == All_Partition )
          {
            if( grid_.comm().size() > 1 )
            {
              // make sure that also ghosts have indices
              markAllUsed<Ghost_Partition>();
            }
          }
    #endif
        }
      }

      //! \copydoc Dune::Fem::ConsecutiveIndexSet::compress */
      bool compress ();

    public:
      ////////////////////////////////////////////////////////
      //  index methods
      //  --index
      ///////////////////////////////////////////////////////
      //! \copydoc Dune::Fem::IndexSet::size */
      template< class Entity >
      IndexType index ( const Entity &entity ) const
      {
        return index< Entity :: codimension >( entity );
      }

      //! \copydoc Dune::Fem::IndexSet::size */
      template< int codim >
      IndexType
      index ( const typename GridPartType::template Codim< codim >::EntityType &entity ) const
      {
        if( codimAvailable( codim ) )
        {
          if( (codim != 0) && ! codimUsed( codim ) )
            setupCodimSet< codim >(std::integral_constant<bool,Dune::Capabilities::hasEntity < GridType, codim > :: v>());

          return codimLeafSet( codim ).index( gridEntity( entity ) );
        }
        else
        {
          DUNE_THROW( NotImplemented, (name() + " does not support indices for codim = ") << codim );
          return -1;
        }
      }

      /* \brief return index for intersection */
      IndexType index ( const IntersectionType &intersection ) const
      {
        enum { codim = intersectionCodimension };
        if( codimAvailable( codim ) )
        {
          // this in only done on first call
          setupIntersections();

          // get corresponding face entity pointer
          FaceType face = getIntersectionFace( intersection );

          return codimLeafSet( codim ).index( gridEntity( face ) );
        }
        else
        {
          DUNE_THROW( NotImplemented, (name() + " does not support indices for intersections, intersectionCodim = ") << codim );
          return -1;
        }
      }

      /* \brief return index for sub entity of given intersection and subEntityNumber */
      IndexType
      subIndex ( const IntersectionType &intersection,
                 int subNumber, unsigned int codim ) const
      {
        DUNE_THROW( NotImplemented, (name() + " does not support subIndices for intersections, intersectionCodim = ") << codim );
        return -1;
      }

      //! \copydoc Dune::Fem::IndexSet::subIndex */
      template< class Entity >
      IndexType subIndex ( const Entity &entity, int subNumber, unsigned int codim ) const
      {
        return subIndex< Entity::codimension >( entity, subNumber, codim );
      }

      //! \copydoc Dune::Fem::IndexSet::subIndex */
      template< int cd >
      IndexType subIndex ( const typename GridPartType::template Codim< cd >::EntityType &entity,
                           int subNumber, unsigned int codim ) const
      {
        assert( (int( codim ) >= cd) && (int( codim ) <= dimension) );
        if( !codimAvailable( codim ) )
          DUNE_THROW( NotImplemented, (name() + " does not support indices for codim = ") << codim );

        if( (codim != 0) && ! codimUsed( codim ) )
        {
          Fem::ForLoop< CallSetUpCodimSet, 1, dimension >::apply( codim, *this );
        }

        const CodimIndexSetType &codimSet = codimLeafSet( codim );
        const IndexType idx = codimSet.subIndex( gridEntity( entity ), subNumber );
        assert( (idx >= 0) && (idx < IndexType( codimSet.size() )) );
        return idx;
      }

      ///////////////////////////////////////////////////////////////////
      //
      //  DoF adjustment methods, needed by AdaptiveDofMapper interface
      //
      ///////////////////////////////////////////////////////////////////

      //! \copydoc Dune::Fem::AdaptiveIndexSet::numberOfHoles */
      IndexType numberOfHoles ( GeometryType type ) const
      {
        const int codim = dimension - type.dim();
        // this index set works only with single geometry types adaptively
        assert( hasSingleGeometryType || geomTypes( codim ).size() == 1 );

        // if mapper is asking for types that are
        // not in the index set then 0 should be returned
        if( geomTypes( codim )[ 0 ] != type )
        {
          return 0;
        }

        return numberOfHoles( codim );
      }

      //! return number of holes of the sets indices
      IndexType numberOfHoles ( const int codim ) const
      {
        if( codimAvailable( codim ) && codimUsed( codim ) )
        {
          return codimLeafSet( codim ).numberOfHoles();
        }
        else
          return 0;
      }

      //! \copydoc Dune::Fem::AdaptiveIndexSet::oldIndex */
      IndexType oldIndex ( IndexType hole, GeometryType type ) const
      {
        const int codim = dimension - type.dim();
        assert( hasSingleGeometryType || geomTypes( codim ).size() == 1 );
        return oldIndex( hole, codim );
      }

      //! return old index, for dof manager only
      IndexType oldIndex (const IndexType hole, const int codim ) const
      {
        if( codimAvailable( codim ) )
        {
          assert( codimUsed( codim ) );
          return codimLeafSet( codim ).oldIndex( hole );
        }
        else
        {
          DUNE_THROW( NotImplemented, (name() + " does not support indices for codim = ") << codim );
          return IndexType(-1);
        }
      }

      //! \copydoc Dune::Fem::AdaptiveIndexSet::newIndex */
      IndexType newIndex ( IndexType hole, GeometryType type ) const
      {
        const int codim = dimension - type.dim();
        assert( hasSingleGeometryType || geomTypes( codim ).size() == 1 );
        return newIndex( hole, codim );
      }

      //! return new index, for dof manager only returns index
      IndexType newIndex (const IndexType hole , const int codim ) const
      {
        if( codimAvailable( codim ) )
        {
          assert( codimUsed( codim ) );
          return codimLeafSet( codim ).newIndex( hole );
        }
        else
        {
          DUNE_THROW( NotImplemented, (name() + " does not support indices for codim = ") << codim );
          return IndexType(-1);
        }
      }

    protected:
      // Note: The following methods forward to Dune::Fem::CodimIndexSet, which
      // expects a Dune::Grid as template argument; all arguments passed to
      // members of Dune::Fem::CodimIndexSet must be compatible with the
      // template grid type. If entities returned by the grid and the grid part
      // respectively differ in type, Dune::Fem::AdaptiveLeafIndexSetBase will
      // call the necessary operations from grid part entities to grid entites.

      // memorise index
      void insertIndex ( const GridElementType &entity );

      // memorise indices for all intersections
      void insertIntersections ( const GridElementType &entity ) const;

      // insert index temporarily
      void insertTemporary ( const GridElementType &entity );

      // set indices to unsed so that they are cleaned on compress
      void removeIndex ( const GridElementType &entity );

      // check whether entity can be inserted or not
      void checkHierarchy ( const GridElementType &entity, bool wasNew );

      // mark indices that are still used (and give new indices to new elements)
      template <PartitionIteratorType pt>
      void markAllUsed ();

      //! clear index set (only for structured grids)
      void clear();

      //! mark all indices of interest
      void setupIndexSet ();

      // give all entities that lie below the old entities new numbers
      // here we need the hierarchic iterator because for example for some
      // grid more the one level of new elements can be created during adaption
      // there for we start to give new number for all elements below the old
      // element
      template <PartitionIteratorType pt>
      void markAllBelowOld ();

      // mark indices that are still used (and give new indices to new elements)
      template< int codim >
      void setupCodimSet (const std::integral_constant<bool,true> &hasEntities) const;
      template< int codim >
      void setupCodimSet (const std::integral_constant<bool,false> &hasEntities) const;

      // mark indices that are still used (and give new indices to new intersections)
      void setupIntersections () const;

      // count elements by iterating over grid and compare
      // entities of given codim with given type
      template< int codim >
      inline IndexType countElements ( GeometryType type, const std::integral_constant<bool,true> &hasEntities ) const;
      template< int codim >
      inline IndexType countElements ( GeometryType type, const std::integral_constant<bool,false> &hasEntities ) const;

    public:
      //! \copydoc Dune::Fem::ConsecutiveIndexSet::write */
      template< class StreamTraits >
      bool write( OutStreamInterface< StreamTraits >& out ) const;

      //! \copydoc Dune::Fem::ConsecutiveIndexSet::read */
      template< class StreamTraits >
      bool read( InStreamInterface< StreamTraits >& in );

    protected:
      FaceType getIntersectionFace( const IntersectionType& intersection ) const
      {
        ElementType inside = intersection.inside();
        return getIntersectionFace( intersection, inside );
      }

      FaceType getIntersectionFace( const IntersectionType& intersection,
                                    const ElementType& inside ) const
      {
        if( ! intersection.conforming() && intersection.neighbor() )
        {
          ElementType outside = intersection.outside();
          // only if outside is more refined then inside
          if( inside.level() < outside.level() )
            return GetFaceEntity :: subEntity( outside, intersection.indexInOutside() );
        }

        // default: get subentity of inside
        return GetFaceEntity :: subEntity( inside, intersection.indexInInside() );
      }
    };

    template< class TraitsImp >
    inline void
    AdaptiveIndexSetBase< TraitsImp >::resizeVectors ()
    {
      codimLeafSet( 0 ).resize();

      // if more than one codimension is supported
      if( numCodimensions > 1 )
      {
        for( int codim = 1; codim < numCodimensions; ++codim )
        {
          if( codimUsed( codim ) )
            codimLeafSet( codim ).resize();
        }
      }
    }


    // --compress
    template< class TraitsImp >
    inline bool
    AdaptiveIndexSetBase< TraitsImp >::compress ()
    {
      // reset list of holes in any case
      for( int codim = 0; codim < numCodimensions; ++codim )
      {
        if( codimUsed( codim ) )
          codimLeafSet( codim ).clearHoles();
      }

      if( compressed_ )
      {
        // if set already compress, do noting for serial runs
        // in parallel runs check sequence number of dof manager
        if( (grid_.comm().size() == 1) || (sequence_ == dofManager_.sequence()) )
          return false;
      }

      // prepare index sets for setup
      for( int codim = 0; codim < numCodimensions; ++codim )
      {
        if( codimUsed( codim ) )
          codimLeafSet( codim ).prepareCompress();
      }

      // mark all indices still needed
      setupIndexSet();

      // true if a least one index is moved
      bool haveToCopy = codimLeafSet( 0 ).compress();
      for( int codim = 1; codim < numCodimensions; ++codim )
      {
        if( codimUsed( codim ) )
          haveToCopy |= codimLeafSet( codim ).compress();
      }

      // now status is compressed
      compressed_ = true;
      // update sequence number
      sequence_ = dofManager_.sequence();

      return haveToCopy;
    }


    template< class TraitsImp >
    inline void
    AdaptiveIndexSetBase< TraitsImp >::insertIndex ( const GridElementType &entity )
    {
#if HAVE_MPI
      // we need special treatment for ghosts
      // ghosts should not be inlcuded in holes list
      if( entity.partitionType() == GhostEntity )
      {
        codimLeafSet( 0 ).insertGhost( entity );
        const bool skipGhosts = (pitype != All_Partition);
        // only for index sets upporting more than one codim
        if( numCodimensions > 1 )
          Fem::ForLoop< InsertGhostSubEntities, 1, dimension >::apply( *this, entity, skipGhosts );
      }
      else
#endif // HAVE_MPI
      {
        codimLeafSet( 0 ).insert( entity );
        // only for index sets supporting more than one codim
        if( numCodimensions > 1 )
          Fem::ForLoop< InsertSubEntities, 1, dimension >::apply( *this, entity );

      }

      assert( codimLeafSet( 0 ).exists( entity ) );

      // insert intersections if this is enabled
      if( intersectionCodimension > 0 )
      {
        insertIntersections( entity );
      }

      // now consecutivity is no longer guaranteed
      compressed_ = false;
    }

    template< class TraitsImp >
    inline void
    AdaptiveIndexSetBase< TraitsImp >::insertIntersections ( const GridElementType &gridElement ) const
    {
      codimLeafSet( intersectionCodimension ).resize();

      const ElementType &element = gridPart_.convert( gridElement );
      const IntersectionIteratorType endiit = gridPart_.iend( element );
      for( IntersectionIteratorType iit = gridPart_.ibegin( element );
           iit != endiit ; ++ iit )
      {
        // get intersection
        const IntersectionType& intersection = *iit ;

        // get correct face pointer
        FaceType face = getIntersectionFace( intersection, element );

        // insert face into index set
        codimLeafSet( intersectionCodimension ).insert( gridEntity( face ) );
      }
    }

    template< class TraitsImp >
    inline void
    AdaptiveIndexSetBase< TraitsImp >::insertTemporary( const GridElementType &entity )
    {
      insertIndex( entity );
      codimLeafSet( 0 ).markForRemoval( entity );
    }

    template< class TraitsImp >
    inline void
    AdaptiveIndexSetBase< TraitsImp >::removeIndex( const GridElementType &entity )
    {
      // remove entities (only mark them as unused)
      codimLeafSet( 0 ).markForRemoval( entity );

      // don't remove higher codim indices (will be done on compression

      // now consecutivity is no longer guaranteed
      compressed_ = false;
    }


    template< class TraitsImp >
    inline void
    AdaptiveIndexSetBase< TraitsImp >
      ::checkHierarchy ( const GridElementType &entity, bool wasNew )
    {
      bool isNew = wasNew ;
      typedef typename GridType::HierarchicIterator HierarchicIterator;

      // for leaf entities, just insert the index
      if( entity.isLeaf() )
      {
        insertIndex( entity );
        return;
      }

      if( isNew )
      {
        // this is a new entity, so insert it,
        // but only temporarily because it's not a leaf entity
        insertTemporary( entity );
      }
      else
      {
        // if we were a leaf entity, all children are new
        isNew = codimLeafSet( 0 ).validIndex( entity );
      }

      // entity has children and we need to go deeper
      const int childLevel = entity.level() + 1;
      const HierarchicIterator end  = entity.hend( childLevel );
      for( HierarchicIterator it = entity.hbegin( childLevel ); it != end; ++it )
        checkHierarchy( *it, isNew );
    }


    template< class TraitsImp >
    template< PartitionIteratorType pt >
    inline void
    AdaptiveIndexSetBase< TraitsImp >::markAllUsed ()
    {
      // make correct size of vectors
      resizeVectors();

      // mark all indices as unused
      for( int codim = 0; codim < numCodimensions; ++codim )
      {
        if( codimUsed( codim ) )
          codimLeafSet( codim ).resetUsed();
      }

      typedef typename GridPartType
        ::template Codim< 0 > :: template Partition< pt > :: IteratorType  Iterator;

      const Iterator end  = gridPart_.template end< 0, pt >();
      for( Iterator it = gridPart_.template begin< 0, pt >(); it != end; ++it )
        insertIndex( gridEntity( *it ) );
    }

    template< class TraitsImp >
    inline void
    AdaptiveIndexSetBase< TraitsImp >::clear()
    {
      // for structured grids clear all information
      // this in only done when setting up grids or after
      // read of parallel data on serial grids
      if( CartesianNonAdaptiveGrid )
      {
        // mark all indices as unused
        for( int codim = 0; codim < numCodimensions; ++codim )
        {
          if( codimUsed( codim ) )
          {
            // clear all information
            codimLeafSet( codim ).clear();
          }
        }
      }
    }

    template< class TraitsImp >
    inline void
    AdaptiveIndexSetBase< TraitsImp >::setupIndexSet ()
    {
      // only done for structured grids
      clear();

#if HAVE_MPI
      // for YaspGrid we need all interior indices first
      if( CartesianNonAdaptiveGrid  &&
          grid_.comm().size() > 1 )
      {
        // we should only get here for YaspGrid
        markAllUsed<Interior_Partition> ();
        if( pitype > Interior_Partition )
          markAllUsed< pitype >();
      }
      else
#endif
      {
        // give all entities that lie on the leaf level new numbers
        markAllUsed< pitype > ();
      }
    }

    template< class TraitsImp >
    template< PartitionIteratorType pt >
    inline void
    AdaptiveIndexSetBase< TraitsImp >::markAllBelowOld ()
    {
      // mark all indices as unused
      for( int codim = 0; codim < numCodimensions; ++codim )
      {
        if( codimUsed( codim ) )
          codimLeafSet( codim ).resetUsed();
      }

      // get macro iterator
      typedef typename GridType::LevelGridView LevelGridView;
      LevelGridView macroView = grid_.levelGridView( 0 );

      typedef typename LevelGridView::template Codim< 0 >::template Partition< pt >::Iterator Iterator;
      const Iterator macroend = macroView.template end< 0, pt >();
      for( Iterator macroit = macroView.template begin< 0, pt >();
           macroit != macroend; ++macroit )
        checkHierarchy( *macroit, false );
    }


    template< class TraitsImp >
    template< int codim >
    inline void
    AdaptiveIndexSetBase< TraitsImp >::setupCodimSet (const std::integral_constant<bool,true>&) const
    {
      // if codim is not available do nothing
      if( ! codimAvailable( codim ) ) return ;

      // create codimLeafSet if not existing
      if( ! codimLeafSet_[ codim ] )
      {
        codimLeafSet_[ codim ].reset( new CodimIndexSetType( grid_, codim ) );
      }

      // resize if necessary
      codimLeafSet( codim ).resize();

      // walk over grid parts entity set and insert entities
      typedef typename GridPartType
        ::template Codim< codim >::template Partition< pitype > :: IteratorType Iterator;

      const Iterator end = gridPart_.template end< codim, pitype >();
      for( Iterator it = gridPart_.template begin< codim, pitype >(); it != end; ++it )
        codimLeafSet( codim ).insert( gridEntity( *it ) );
    }

    template< class TraitsImp >
    template< int codim >
    inline void
    AdaptiveIndexSetBase< TraitsImp >::setupCodimSet (const std::integral_constant<bool,false>&) const
    {
      // if codim is not available do nothing
      if( ! codimAvailable( codim ) ) return ;

      // create codimLeafSet if not existing
      if( ! codimLeafSet_[ codim ] )
      {
        codimLeafSet_[ codim ].reset( new CodimIndexSetType( grid_, codim ) );
      }

      // resize if necessary
      codimLeafSet( codim ).resize();

      typedef typename GridPartType
        ::template Codim< 0 >::template Partition< pitype > :: IteratorType Iterator;

      const Iterator end = gridPart_.template end< 0, pitype >();
      for( Iterator it = gridPart_.template begin< 0, pitype >(); it != end; ++it )
      {
        const ElementType& element = *it ;
        const GridElementType &gridElement = gridEntity( element );
        InsertSubEntities< codim >::apply( *this, gridElement );
      }
    }


    template< class TraitsImp >
    inline void
    AdaptiveIndexSetBase< TraitsImp >::setupIntersections() const
    {
      // if intersectionCodimension < 0 then this feature is disabled
      if( intersectionCodimension < 0 ) return ;

      // do nothing if insections are already available
      if( codimUsed( intersectionCodimension ) ) return ;

      // resize if necessary
      codimLeafSet( intersectionCodimension ).resize();

      // walk over grid parts entity set and insert entities
      typedef typename GridPartType
        ::template Codim< 0 >::template Partition< pitype > :: IteratorType Iterator;

      const Iterator end = gridPart_.template end< 0, pitype >();
      for( Iterator it = gridPart_.template begin< 0, pitype >(); it != end; ++it )
      {
        // insert all intersections of this entity
        insertIntersections( gridEntity( *it ) );
      }
    }

    template< class TraitsImp >
    template< int codim >
    inline typename AdaptiveIndexSetBase< TraitsImp >::IndexType
    AdaptiveIndexSetBase< TraitsImp >::countElements ( GeometryType type, const std::integral_constant<bool,true>& ) const
    {
      typedef typename GridPartType
        ::template Codim< codim > :: template Partition< pitype > :: IteratorType Iterator;

      const Iterator begin = gridPart_.template begin< codim, pitype >();
      const Iterator end = gridPart_.template end< codim, pitype >();
      IndexType count = 0;
      for( Iterator it = begin; it != end; ++it )
      {
        if( it->type() == type )
        {
          ++count;
        }
      }
      return count;
    }

    template< class TraitsImp >
    template< int codim >
    inline typename AdaptiveIndexSetBase< TraitsImp >::IndexType
    AdaptiveIndexSetBase< TraitsImp >::countElements ( GeometryType type, const std::integral_constant<bool,false>& ) const
    {
      if( ! codimLeafSet_[ codim ] )
        return 0;

      // make sure codimension is enabled
      assert( codimAvailable( codim ) );

      // resize if necessary
      codimLeafSet( codim ).resize();

      typedef typename GridPartType
        ::template Codim< 0 >::template Partition< pitype > :: IteratorType Iterator;

      typedef typename GridPartType::ctype ctype;

      const Iterator end = gridPart_.template end< 0, pitype >();
      IndexType count = 0;
      for( Iterator it = gridPart_.template begin< 0, pitype >(); it != end; ++it )
      {
        const ElementType& element = *it ;
        const GridElementType &gridElement = gridEntity( element );
        const int subEntities = gridElement.subEntities( codim );
        for (int i=0; i < subEntities; ++i)
        {
          if (! codimLeafSet( codim ).exists( gridElement, i) )
          {
            codimLeafSet( codim ).insertSubEntity( gridElement,i );
            if ( Dune::ReferenceElements< ctype, dimension >::
               general( gridElement.type() ).type( i, codim ) == type )
            {
              ++count;
            }
          }
        }
      }

      return count;
    }


    template< class TraitsImp >
    template< class StreamTraits >
    inline bool AdaptiveIndexSetBase< TraitsImp >
      ::write ( OutStreamInterface< StreamTraits >& out ) const
    {
      // write name for indentification
      const std::string myname( name() );
      out << myname;

      // write number of codimensions
      out << numCodimensions ;

      // write whether codim is used
      for( int i = 0; i < numCodimensions; ++i )
        out << codimUsed( i );

      // write all sets
      for( int i = 0; i < numCodimensions; ++i )
      {
        if( codimUsed( i ) )
          codimLeafSet( i ).write( out );
      }

      // if we got until here writing was sucessful
      return true;
    }


    template< class TraitsImp >
    template< class StreamTraits >
    inline bool AdaptiveIndexSetBase< TraitsImp >
      ::read ( InStreamInterface< StreamTraits > &in )
    {
      {
        // read name and check compatibility
        std::string storedName;
        in >> storedName;

        std::string myname( name() );

        if( myname != storedName )
        {
          size_t length = std::min( myname.size(), storedName.size() );
          // only print the first character of whatever storedName is
          std::string found = storedName.substr(0, length-1 );
          DUNE_THROW( InvalidStateException,
                      "AdaptiveIndexSetBase::read: got  " << found
                      << " (expected " << myname << ")." );
        }
      }

      // read number of codimensions
      int numCodim;
      in >> numCodim;

      // make sure everything is correct
      // assert( numCodim == numCodimensions );
      if( numCodim != numCodimensions )
        DUNE_THROW(InvalidStateException,"AdaptiveIndexSetBase::read: got wrong number of codimensions" << numCodim << " instead of " << numCodimensions);

      // read codim used
      for( int i = 0; i < numCodimensions; ++i )
      {
        bool codimInUse = false ;
        in >> codimInUse;
        if( codimInUse && ! codimLeafSet_[ i ] )
        {
          codimLeafSet_[ i ].reset( new CodimIndexSetType( grid_, (i == intersectionCodimension ) ? 1 : i  ) );
        }
      }

      for( int i = 0; i < numCodimensions; ++i )
      {
        if( codimUsed( i ) )
          codimLeafSet( i ).read( in );
      }

      // in parallel runs we have to compress here
      if( grid_.comm().size() > 1 )
        compressed_ = false;

      // if we got until here reading was sucessful
      return true;
    }



    /////////////////////////////////////////////////////////////////////////
    //
    //  --AdaptiveLeafIndexSet
    //
    /////////////////////////////////////////////////////////////////////////

    template< class GridPartImp >
    struct AdaptiveLeafIndexSetTraits
    : public AdaptiveIndexSetBaseTraits< GridPartImp, AdaptiveLeafIndexSet< GridPartImp > >
    {
      // number of codimensions
      enum { numCodimensions = GridPartImp :: dimension + 1 };
      // first comdimension that is supported (not yet supported)
      enum { startingCodimension = 0 };
      // intersection codimensions (this is usually dimension + 1 )
      enum { intersectionCodimension = -1 };
    };

    /** \class AdaptiveLeafIndexSet
     *  \brief consecutive, persistent index set for the leaf level based on the
     *         grid's hierarchy index set
     *
     *  This index set generates a consecutive leaf index out of the unique global
     *  index of each entity. It can be used instead of the default grid index sets
     *  and can be generated for each grid implementation.
     *
     *  \note This index sets supports all indices for all codimensions of the grid.
     *
     */
    template < class GridPartImp >
    class AdaptiveLeafIndexSet
    : public AdaptiveIndexSetBase< AdaptiveLeafIndexSetTraits< GridPartImp > >
    {
      typedef AdaptiveIndexSetBase< AdaptiveLeafIndexSetTraits< GridPartImp > > BaseType;
      typedef AdaptiveLeafIndexSetTraits< GridPartImp > Traits;

    public:
      typedef typename BaseType :: GridType GridType;
      typedef typename BaseType :: GridPartType GridPartType;

      //! Constructor
      AdaptiveLeafIndexSet (const GridType* grid)
        : BaseType(grid)
      {}

      //! Constructor
      AdaptiveLeafIndexSet (const GridPartType& gridPart)
        : BaseType(gridPart)
      {}

      //! return name of index set
      virtual std::string name () const
      {
        return "AdaptiveLeafIndexSet";
      }

      bool compress ()
      {
        const bool compressed = BaseType::compress();

#ifndef NDEBUG
        if( this->grid_.comm().size() == 1 )
        {
          for( int codim = Traits::startingCodimension; codim < Traits::numCodimensions; ++codim )
            assert( this->codimUsed( codim ) ?
                size_t(this->size( codim )) == size_t(this->grid_.size( codim )) : true );
        }
#endif // #ifndef NDEBUG

        return compressed;
      }
    };


    /////////////////////////////////////////////////////////////////////////
    //
    //  --IntersectionAdaptiveLeafIndexSet
    //
    /////////////////////////////////////////////////////////////////////////

    template< class GridPartImp >
    struct IntersectionAdaptiveLeafIndexSetTraits
    : public AdaptiveIndexSetBaseTraits< GridPartImp, IntersectionAdaptiveLeafIndexSet< GridPartImp > >
    {
      // number of codimensions
      enum { numCodimensions = GridPartImp :: dimension + 2 };
      // intersection codimensions (this is usually dimension + 1 )
      enum { intersectionCodimension = numCodimensions - 1 };
      // first comdimension that is supported (not yet supported)
      enum { startingCodimension = 0 };
    };

    /** \class AdaptiveLeafIndexSet
     *  \brief consecutive, persistent index set for the leaf level based on the
     *         grid's hierarchy index set
     *
     *  This index set generates a consecutive leaf index out of the unique global
     *  index of each entity. It can be used instead of the default grid index sets
     *  and can be generated for each grid implementation.
     *
     *  \note This index sets supports all indices for all codimensions of the grid.
     *
     */
    template < class GridPartImp >
    class IntersectionAdaptiveLeafIndexSet
    : public AdaptiveIndexSetBase< IntersectionAdaptiveLeafIndexSetTraits< GridPartImp > >
    {
      typedef AdaptiveIndexSetBase< IntersectionAdaptiveLeafIndexSetTraits< GridPartImp > > BaseType;
      typedef IntersectionAdaptiveLeafIndexSetTraits< GridPartImp > Traits;

    public:
      typedef typename BaseType :: GridType     GridType;
      typedef typename BaseType :: GridPartType GridPartType;
      //! Constructor
      IntersectionAdaptiveLeafIndexSet (const GridType* grid)
        : BaseType(grid)
      {}

      //! Constructor
      IntersectionAdaptiveLeafIndexSet (const GridPartType& gridPart)
        : BaseType(gridPart)
      {}

      //! return name of index set
      virtual std::string name () const
      {
        return "IntersectionAdaptiveLeafIndexSet";
      }

      bool compress ()
      {
        const bool compressed = BaseType::compress();

#ifndef NDEBUG
        if( this->grid_.comm().size() == 1 )
        {
          for( int codim = Traits::startingCodimension; codim < Traits::numCodimensions; ++codim )
            if( codim != Traits::intersectionCodimension )
              assert( size_t(this->size( codim )) == size_t(this->grid_.size( codim )) );
        }
#endif // #ifndef NDEBUG

        return compressed;
      }
    };

    /////////////////////////////////////////////////////////////////////////
    //
    //  --DGAdaptiveLeafIndexSet
    //
    /////////////////////////////////////////////////////////////////////////

    template< class GridPartImp >
    struct DGAdaptiveLeafIndexSetTraits
    : public AdaptiveIndexSetBaseTraits< GridPartImp, DGAdaptiveLeafIndexSet< GridPartImp > >
    {
      // this index set only supports one codimension, codim zero
      enum { numCodimensions = 1 };
      // first comdimension that is supported (not yet supported)
      enum { startingCodimension = 0 };
      // intersection codimensions (this is usually dimension + 1 )
      enum { intersectionCodimension = -1 };
    };

    /** \class DGAdaptiveLeafIndexSet
     *  \brief consecutive, persistent index set for the leaf level based on the
     *         grid's hierarchy index set
     *
     *  This index set generates a consecutive leaf index out of the unique global
     *  index of each codimension 0 entity.
     *
     *  \note This index sets supports only indices for codimensions 0 entities of the grid.
     *
     */
    template < class GridPartImp >
    class DGAdaptiveLeafIndexSet
    : public AdaptiveIndexSetBase< DGAdaptiveLeafIndexSetTraits< GridPartImp > >
    {
      typedef AdaptiveIndexSetBase< DGAdaptiveLeafIndexSetTraits< GridPartImp > > BaseType;
      typedef DGAdaptiveLeafIndexSetTraits< GridPartImp > Traits;

    public:
      typedef typename BaseType :: GridType     GridType;
      typedef typename BaseType :: GridPartType GridPartType;
      //! Constructor
      DGAdaptiveLeafIndexSet (const GridType* grid)
        : BaseType(grid)
      {}

      DGAdaptiveLeafIndexSet (const GridPartType& gridPart)
        : BaseType(gridPart)
      {}

      //! return name of index set
      virtual std::string name () const
      {
        return "DGAdaptiveLeafIndexSet";
      }

      bool compress ()
      {
        const bool compressed = BaseType::compress();

#ifndef NDEBUG
        if( this->grid_.comm().size() == 1 )
          assert( size_t(this->size( 0 )) == size_t(this->grid_.size( 0 )) );
#endif // #ifndef NDEBUG

        return compressed;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ADAPTIVELEAFINDEXSET_HH
