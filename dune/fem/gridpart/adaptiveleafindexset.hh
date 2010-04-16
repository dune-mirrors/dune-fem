#ifndef DUNE_ADAPTIVELEAFINDEXSET_HH
#define DUNE_ADAPTIVELEAFINDEXSET_HH

//- local includes 
#include <dune/common/forloop.hh>
#include <dune/fem/gridpart/dunefemindexsets.hh>
#include <dune/fem/gridpart/codimindexset.hh>

namespace Dune
{
  /** \class AdaptiveLeafIndexSetBase
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
  class AdaptiveLeafIndexSetBase
  : public ConsecutivePersistentIndexSet< 
            typename TraitsImp :: GridType, 
            AdaptiveLeafIndexSetBase< TraitsImp > 
            >
  {
  protected:  
    typedef typename TraitsImp :: GridType GridType;
    typedef AdaptiveLeafIndexSetBase< TraitsImp > ThisType;
    typedef ConsecutivePersistentIndexSet< GridType, ThisType > BaseType;

    friend class Conversion< ThisType, EmptyIndexSet >;

  public:
    //! dimension of the grid 
    static const int dimension = GridType::dimension;

    //! number of supported codimensions 
    static const int numCodimensions = TraitsImp :: numCodimensions ;

    //! type of index 
    typedef typename BaseType :: IndexType IndexType;

    //! type of codimension 0 Entity 
    typedef typename GridType::template Codim< 0 >::Entity ElementType;

  private:
    typedef CodimIndexSet CodimIndexSetType; 

    typedef HierarchicIndexSetSelector< GridType > SelectorType;
    typedef typename SelectorType :: HierarchicIndexSet HIndexSetType;

    template< int codim >
    struct CountElements
    {
      static void apply ( const ThisType &indexSet, const GeometryType &type, int &count )
      {
        if( type.dim() == dimension - codim )
          count = indexSet.template countElements< codim >( type );
      }
    };

    template< int codim >
    struct InsertSubEntities
    {
      static void apply ( ThisType &indexSet, const ElementType &entity )
      {
        // if codimension is not available return 
        if( ! indexSet.codimAvailable( codim ) ) return ;

        // if codimension is not used return 
        if( !indexSet.codimUsed_[ codim ] ) return;
        
        const HIndexSetType &hIndexSet = indexSet.hIndexSet_;
        CodimIndexSetType &codimSet = indexSet.codimLeafSet_[ codim ];

        for( int i = 0; i < entity.template count< codim >(); ++i )
        {
          // codimSet.insert( hIndexSet.template subIndex< codim >( entity, i ) );
          codimSet.insert( hIndexSet.subIndex( entity, i, codim ) );
        }
      }
    };

    template< int codim >
    struct InsertGhostSubEntities
    {
      static void apply ( ThisType &indexSet, const ElementType &entity ,
                          const bool skipGhosts )
      {
        // if codimension is not available return 
        if( ! indexSet.codimAvailable( codim ) ) return ;

        // if codimension is not used return 
        if( !indexSet.codimUsed_[ codim ] ) return;

        const HIndexSetType &hIndexSet = indexSet.hIndexSet_;
        CodimIndexSetType &codimSet = indexSet.codimLeafSet_[ codim ];

        for( int i = 0; i < entity.template count< codim >(); ++i )
        {
          typedef typename GridType::template Codim< codim >::EntityPointer EntityPointer;
          typedef typename GridType::template Codim< codim >::Entity Entity;

          EntityPointer ptr = entity.template subEntity< codim >( i );
          const Entity &subentity = *ptr;

          if( !skipGhosts || (entity.partitionType() != GhostEntity) )
            codimSet.insertGhost( hIndexSet.index( subentity ) );
        }
      }
    };

    template< int codim >
    struct CallSetUpCodimSet
    {
      static void apply ( const int cd, const ThisType &indexSet )
      {
        // if codimension is not available return 
        if( ! indexSet.codimAvailable( codim ) ) return ;

        if( cd == codim )
          indexSet.template setUpCodimSet< codim >();
      }
    };

    //! is true if grid is structured grid 
    enum { StructuredGrid = ! Capabilities::IsUnstructured<GridType>::v };

    // my type, to be revised 
    enum { myType = ( numCodimensions == 1 ) ? ( (StructuredGrid) ? -1 : 665 ) : 6 };
    enum { myVersionTag = -665 };

    //! default partition iterator type 
    static const PartitionIteratorType pitype = TraitsImp :: pitype;

    // reference to HierarchicIndexSet 
    const HIndexSetType &hIndexSet_;
    // Codimension leaf index sets 
    mutable CodimIndexSetType codimLeafSet_[ numCodimensions ];
    // flag for codim is in use or not 
    mutable bool codimUsed_ [ numCodimensions ];
    
    // actual sequence number 
    int sequence_;

    //! flag is tru if set is in compressed status
    mutable bool compressed_;

  protected:
    using BaseType::grid_;
    using BaseType::dofManager_;

    // return true if codim is supported 
    const bool codimAvailable( const int codim ) const 
    {
      return codim < numCodimensions ;
    }

  public:
    //! type traits of this class (see defaultindexsets.hh)
    typedef DefaultLeafIteratorTypes<GridType> Traits; 

    //! Constructor
    AdaptiveLeafIndexSetBase (const GridType & grid) 
      : BaseType(grid) 
      , hIndexSet_( SelectorType::hierarchicIndexSet(grid) ) 
      , sequence_( this->dofManager_.sequence() )
      , compressed_(true) // at start the set is compressed 
    {
      // codim 0 is used by default
      codimUsed_[ 0 ] = true;

      // all higher codims are not used by default
      for(int codim = 1; codim < numCodimensions; ++codim ) codimUsed_[ codim ] = false ;
      
      // set the codim of each Codim Set. 
      for(int codim = 0; codim < numCodimensions; ++codim) codimLeafSet_[ codim ].setCodim( codim );

      // build index set 
      setupIndexSet();
    }

    //! Destructor
    virtual ~AdaptiveLeafIndexSetBase ()
    {}

    //! return type of index set, for GrapeDataIO
    int type () const
    {
      return myType;
    }

    //! return name of index set 
    virtual std::string name () const
    {
      return "AdaptiveLeafIndexSetBase";
    }

    //****************************************************************
    //
    //  INTERFACE METHODS for DUNE INDEX SETS 
    //
    //****************************************************************
    //! return size of grid entities per level and codim 
    IndexType size ( GeometryType type ) const
    {
      const int codim = dimension - type.dim();

      // use size of codim index set if possible 
      if( codimAvailable( codim ) )
      {
        if( codimUsed_[ codim ] )
          return codimLeafSet_[ codim ].size();
      }

      assert( hIndexSet_.geomTypes( codim ).size() == 1 );
      int count = 0;
      ForLoop< CountElements, 0, dimension > :: apply( *this, type, count );
      return count;
    }
    
    //! return size of grid entities of given codim 
    IndexType size ( int codim ) const
    {
      assert( hIndexSet_.geomTypes(codim).size() == 1 ); 
      return size(hIndexSet_.geomTypes(codim)[0]);
    }
    
    //! returns vector with geometry tpyes this index set has indices for
    const std::vector <GeometryType> & geomTypes (int codim) const 
    {
      return hIndexSet_.geomTypes(codim);
    }
   
    //! \brief returns true if entity is contained in index set 
    template <class EntityType>
    bool contains (const EntityType & en) const
    {
      enum { codim = EntityType::codimension };
      if( codimAvailable( codim ) )
      {
        assert( codimUsed_[codim] );
        return codimLeafSet_[codim].exists( hIndexSet_.index( en ) ); 
      }
      else 
        return false;
    }

    //****************************************************************
    //
    //  METHODS for Adaptation with DofManger
    //
    //****************************************************************

    //! insert entity into set 
    void insertEntity( const ElementType &entity )  
    {
      // here we have to add the support of higher codims 
      resizeVectors();
      insertIndex( entity );
    }

    //! Unregister entity from set 
    void removeEntity( const ElementType &entity )
    {
      removeIndex( entity );
    }

    //! reallocate the vector for new size
    void resizeVectors ();

    //! if grid has changed, resize index vectors, and create 
    //! indices for new entities, new entities are entities that 
    //! lie below the old entities 
    void resize ()
    {
      resizeVectors();

  #if HAVE_MPI
      if( StructuredGrid  &&
          this->grid_.comm().size() > 1 )
      {
        codimLeafSet_.clear();

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
          if( this->grid_.comm().size() > 1 )
          {
            // make sure that also ghosts have indices 
            markAllUsed<Ghost_Partition>();
          }
        }
  #endif
      }
    }

    //! make to index numbers consecutive 
    //! return true, if at least one hole existed  
    bool compress ();

  public:
    template< class Entity >
    IndexType index ( const Entity &entity ) const
    {
      return index< Entity :: codimension >( entity );
    }

    template< int codim >
    IndexType
    index ( const typename GridType :: template Codim< codim > :: Entity &entity ) const
    {
      if( codimAvailable( codim ) ) 
      {
        if( (codim != 0) && ! codimUsed_[ codim ] )
          setUpCodimSet< codim >();

        //assertCodimSetSize( codim );
        const CodimIndexSetType &codimSet = codimLeafSet_[ codim ];
        const int hIdx = hIndexSet_.template index( entity );
        const int idx = codimSet.index( hIdx );
        assert( (idx >= 0) && (idx < codimSet.size()) );
        return idx;
      }
      else 
      {
        assert( false );
        return -1;
      }
    }

    template< int codim >
    IndexType DUNE_DEPRECATED
    subIndex ( const typename GridType :: template Codim< 0 > :: Entity &entity,
               int subNumber ) const
    {
      if( codimAvailable( codim ) ) 
      {
        if( (codim != 0) && !codimUsed_[ codim ] )
          setUpCodimSet< codim >();

        //assertCodimSetSize( codim );
        const CodimIndexSetType &codimSet = codimLeafSet_[ codim ];
        const int hIdx = hIndexSet_.template subIndex< codim >( entity, subNumber );
        const int idx = codimSet.index( hIdx );
        assert( (idx >= 0) && (idx < codimSet.size()) );
        return idx;
      }
      else 
      {
        assert( false );
        return -1;
      }
    }

    IndexType
    subIndex ( const typename GridType::template Codim< 0 >::Entity &entity,
               int subNumber, unsigned int codim ) const
    {
      if( codimAvailable( codim ) ) 
      {
        if( (codim != 0) && !codimUsed_[ codim ] )
          ForLoop< CallSetUpCodimSet, 0, dimension >::apply( codim, *this );
        
        const CodimIndexSetType &codimSet = codimLeafSet_[ codim ];
        const int hIdx = hIndexSet_.subIndex( entity, subNumber, codim );
        const int idx = codimSet.index( hIdx );
        assert( (idx >= 0) && (idx < codimSet.size()) );
        return idx;
      }
      else 
      {
        assert( false );
        return -1;
      }
    }

    //! return number of holes of the sets indices 
    int numberOfHoles ( const int codim ) const
    {
      if( codimAvailable( codim ) ) 
      {
        assert( codimUsed_[codim] );
        return codimLeafSet_[codim].numberOfHoles(); 
      }
      else 
        return 0;
    }

    //! return old index, for dof manager only 
    int oldIndex (const int hole, const int codim ) const
    {
      if( codimAvailable( codim ) ) 
      {
        assert( codimUsed_[codim] );
        return codimLeafSet_[codim].oldIndex( hole ); 
      }
      else 
      {
        assert( false );
        return -1; 
      }
    }

    //! return new index, for dof manager only returns index 
    int newIndex (const int hole , const int codim ) const
    {
      if( codimAvailable( codim ) ) 
      {
        assert( codimUsed_[codim] );
        return codimLeafSet_[codim].newIndex( hole ); 
      }
      else 
      {
        assert( false );
        return -1; 
      }
    }

  protected:
    // memorise index 
    void insertIndex ( const ElementType &entity );

    // set indices to unsed so that they are cleaned on compress
    void removeIndex ( const ElementType &entity );

    // check whether entity can be inserted or not 
    void checkHierarchy ( const ElementType &entity, const bool wasNew );

    // mark indices that are still used (and give new indices to new elements)
    template <PartitionIteratorType pt>
    void markAllUsed (); 

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
    void setUpCodimSet () const;

    // count elements by iterating over grid and compare 
    // entities of given codim with given type 
    template< int codim >
    inline int countElements ( GeometryType type ) const;
    
  public:
    //! write indexset to stream  
    template< class StreamTraits > 
    bool write( OutStreamInterface< StreamTraits >& out ) const;

    //! read indexset from stream  
    template< class StreamTraits > 
    bool read( InStreamInterface< StreamTraits >& in );

    //! write indexset to xdr file 
    bool write_xdr( const std::string &filename ) const ;
    //! write indexset to xdr file 
    bool write_xdr( const std::string &filename, int timestep ) const ;

    //! read index set from given xdr file 
    bool read_xdr( const std::string &filename, int timestep );
    //! read index set from given xdr file 
    bool read_xdr( const std::string &filename );
  };

  template< class TraitsImp >
  inline void
  AdaptiveLeafIndexSetBase< TraitsImp >::resizeVectors ()
  {
    codimLeafSet_[ 0 ].resize( hIndexSet_.size( 0 ) );

    // if more than one codimension is supported 
    if( numCodimensions > 1 ) 
    {
      for( int codim = 1; codim < numCodimensions; ++codim )
      {
        if( codimUsed_[ codim ] )
          codimLeafSet_[ codim ].resize( hIndexSet_.size( codim ) );
      }
    }
  }


  // --compress
  template< class TraitsImp >
  inline bool
  AdaptiveLeafIndexSetBase< TraitsImp >::compress ()
  {
    // reset list of holes in any case
    for( int codim = 0; codim < numCodimensions; ++codim )
      codimLeafSet_[ codim ].clearHoles();

    if( compressed_ )
    {
      // if set already compress, do noting for serial runs
      // in parallel runs check sequence number of dof manager 
      if( (grid_.comm().size() == 1) || (sequence_ == dofManager_.sequence()) )
        return false;
    }

    // mark all indices still needed 
    setupIndexSet();

    // true if a least one index is moved
    bool haveToCopy = codimLeafSet_[ 0 ].compress();
    for( int codim = 1; codim < numCodimensions; ++codim )
    {
      if( codimUsed_[ codim ] ) 
        haveToCopy |= codimLeafSet_[ codim ].compress();
    }

    // now status is compressed
    compressed_ = true;
    // update sequence number
    sequence_ = dofManager_.sequence();

    return haveToCopy;
  }


  template< class TraitsImp >
  inline void
  AdaptiveLeafIndexSetBase< TraitsImp >::insertIndex ( const ElementType &entity )
  {
    const int index = hIndexSet_.index( entity );

#if HAVE_MPI 
    // we need special treatment for ghosts 
    // ghosts should not be inlcuded in holes list 
    if( entity.partitionType() == GhostEntity )
    {
      codimLeafSet_[ 0 ].insertGhost( index );
      const bool skipGhosts = (pitype != All_Partition);
      // only for index sets upporting more than one codim 
      if( numCodimensions > 1 )
        ForLoop< InsertGhostSubEntities, 1, dimension >::apply( *this, entity, skipGhosts );
    }
    else 
#endif // HAVE_MPI
    {
      codimLeafSet_[ 0 ].insert( index );
      // only for index sets upporting more than one codim 
      if( numCodimensions > 1 )
        ForLoop< InsertSubEntities, 1, dimension >::apply( *this, entity );
    }

    assert( codimLeafSet_[ 0 ].exists( index ) );

    // now consecutivity is no longer guaranteed
    compressed_ = false;
  }


  template< class TraitsImp >
  inline void
  AdaptiveLeafIndexSetBase< TraitsImp >::removeIndex( const ElementType &entity )
  {
    // remove entities (only mark them as unused)
    codimLeafSet_[ 0 ].remove( hIndexSet_.index( entity ) );

    // don't remove higher codim indices (will be done on compression

    // now consecutivity is no longer guaranteed
    compressed_ = false;
  }


  template< class TraitsImp >
  inline void
  AdaptiveLeafIndexSetBase< TraitsImp >
    ::checkHierarchy ( const ElementType &entity, bool isNew )
  {
    typedef typename ElementType::HierarchicIterator HierarchicIterator;

    // for leaf entites, just insert the index
    if( entity.isLeaf() )
    {
      insertIndex( entity );
      return;
    }

    if( isNew )
    {
      // this is a new entity, so insert it
      insertIndex( entity );
      // but it's not a leaf entity, so mark it unused
      removeIndex( entity );
    }
    else
    {
      // if we were a leaf entity, all children are new
      isNew = codimLeafSet_[ 0 ].validIndex( hIndexSet_.index ( entity ) );
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
  AdaptiveLeafIndexSetBase< TraitsImp >::markAllUsed ()
  {
    // make correct size of vectors 
    resizeVectors();

    // mark all indices as unused
    for( int codim = 0; codim < numCodimensions; ++codim )
    {
      if( codimUsed_[ codim ] )
        codimLeafSet_[ codim ].set2Unused();
    }
    
    typedef typename GridType
      ::template Codim< 0 >::template Partition< pt >::LeafIterator
      Iterator;

    const Iterator end  = grid_.template leafend< 0, pt >();
    for( Iterator it = grid_.template leafbegin< 0, pt >(); it != end; ++it )
      insertIndex( *it );
  }

  template< class TraitsImp >
  inline void
  AdaptiveLeafIndexSetBase< TraitsImp >::setupIndexSet ()
  {
    // for structured grids clear all information 
    // this in only done when setting up grids or after 
    // read of parallel data on serial grids 
    if( StructuredGrid )
    {
      // mark all indices as unused
      for( int codim = 0; codim < numCodimensions; ++codim )
      {
        if( codimUsed_[ codim ] )
        {
          // clear all information 
          codimLeafSet_[ codim ].clear();
        }
      }
    }

#if HAVE_MPI
    // for YaspGrid we need all interior indices first 
    // so we can use SGrid for the visualization :(
    if( StructuredGrid  &&
        this->grid_.comm().size() > 1 )
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
  AdaptiveLeafIndexSetBase< TraitsImp >::markAllBelowOld ()
  {
    // mark all indices as unused
    for( int codim = 0; codim < numCodimensions; ++codim )
    {
      if( codimUsed_[ codim ] )
        codimLeafSet_[ codim ].set2Unused(); 
    }
    
    // get macro iterator 
    typedef typename GridType
      ::template Codim< 0 >::template Partition< pt >::LevelIterator
      Iterator;

    const Iterator macroend = grid_.template lend< 0, pt >( 0 );
    for( Iterator macroit = grid_.template lbegin< 0, pt >( 0 ); 
         macroit != macroend; ++macroit )
      checkHierarchy( *macroit, false );
  }


  template< class TraitsImp >
  template< int codim >
  inline void
  AdaptiveLeafIndexSetBase< TraitsImp >::setUpCodimSet () const
  {
    // if codim is not available do nothing 
    if( ! codimAvailable( codim ) ) return ;

    // resize if necessary 
    codimLeafSet_[ codim ].resize( hIndexSet_.size( codim ) );
    
    // walk over leaf level on locate all needed entities  
    typedef typename GridType
      ::template Codim< codim >::template Partition< pitype >::LeafIterator
      Iterator;

    const Iterator end = grid_.template leafend< codim, pitype >();
    for( Iterator it = grid_.template leafbegin< codim, pitype >(); it != end; ++it )
      codimLeafSet_[ codim ].insert( hIndexSet_.index ( *it ) );

    // mark codimension as used
    codimUsed_[ codim ] = true;
  }


  template< class TraitsImp >
  template< int codim >
  inline int
  AdaptiveLeafIndexSetBase< TraitsImp >::countElements ( GeometryType type ) const
  {
    typedef typename GridType
      ::template Codim< codim >::template Partition< pitype >::LeafIterator
      Iterator;

    const GridType &grid = this->grid_;

    int count = 0;
    const Iterator begin = grid.template leafbegin< codim, pitype >();
    const Iterator end = grid.template leafend< codim, pitype >();
    for( Iterator it = begin; it != end; ++it )
    {
      if( it->type() == type )
        ++count;
    }
    return count; 
  }


  template< class TraitsImp >
  template< class StreamTraits > 
  inline bool AdaptiveLeafIndexSetBase< TraitsImp >
    ::write ( OutStreamInterface< StreamTraits >& out ) const
  {
    // write new verion tag 
    int newVersion = myVersionTag;
    out << newVersion ;

    // write my type
    int typeVar = type();
    out << typeVar;

    // write whether codim is used 
    for( int i = 0; i < numCodimensions; ++i )
      out << codimUsed_[ i ];

    // write all sets 
    for( int i = 0; i < numCodimensions; ++i )
      codimLeafSet_[ i ].write( out );
    
    // if we got until here writing was sucessful
    return true;
  }

  template< class TraitsImp >
  inline bool AdaptiveLeafIndexSetBase< TraitsImp >
    ::write_xdr ( const std::string &filename ) const
  {
#if DUNE_FEM_COMPATIBILITY
    // create write stream
    XDRWriteStream xdr( filename );
    bool success = true;

    // write new verion tag 
    int newVerion = myVersionTag;
    success &= xdr.inout( newVerion );

    // write my type
    int typeVar = type();
    success &= xdr.inout( typeVar );

    // write whether codim is used 
    for( int i = 0; i < numCodimensions; ++i )
      success &= xdr.inout( codimUsed_[ i ] );

    // write all sets 
    for( int i = 0; i < numCodimensions; ++i )
      success &= codimLeafSet_[ i ].processXdr( xdr );
    
    return success;
#else 
    // new version using streams 
    try
    {
      XDRFileOutStream out( filename );
      return write( out );
    }
    catch( Exception e )
    {
      return false;
    }
#endif
  }


  template< class TraitsImp >
  inline bool AdaptiveLeafIndexSetBase< TraitsImp >
    ::write_xdr ( const std::string &filename, int timestep ) const
  {
    const char *path = "";
    std::string fnstr = genFilename( path, filename, timestep );
    return write_xdr( fnstr );
  }


  template< class TraitsImp >
  template< class StreamTraits > 
  inline bool AdaptiveLeafIndexSetBase< TraitsImp >
    ::read ( InStreamInterface< StreamTraits > &in ) 
  {
    // check new version tag
    int newVersionTag = myVersionTag;
    in >> newVersionTag;
    const bool newVersion = (newVersionTag == myVersionTag);

    // if new version the read type, otherwise newVersionTag is the type info
    int typeVar = (newVersion ? type() : newVersionTag);
    if( newVersion )
      in >> typeVar;

    // index set type check
    if( (typeVar != 2) && (typeVar != type()) )
    {
      DUNE_THROW( InvalidStateException,
                  "AdaptiveLeafIndexSetBase::read: wrong type " << typeVar
                  << " given (expected " << type() << ")." );
    }

    if( newVersion )
    {
      // read codim used 
      for( int i = 0; i < numCodimensions; ++i )
        in >> codimUsed_[ i ];
    }
    else 
    {
      // it depends on the type whether higher codims were stored
      for( int i = 0; i < numCodimensions; ++i )
        codimUsed_[ i ] = (typeVar == type());
    }

    if( typeVar == type() )
    {
      for( int i = 0; i < numCodimensions; ++i )
      {
        if( codimUsed_[ i ] )
          codimLeafSet_[ i ].read( in );
      }
    }
    else
      codimLeafSet_[ 0 ].read( in );

    // in parallel runs we have to compress here
    if( grid_.comm().size() > 1 )
      compressed_ = false;
    
    // if we got until here reading was sucessful 
    return true;
  }

  template< class TraitsImp >
  inline bool AdaptiveLeafIndexSetBase< TraitsImp >
    ::read_xdr ( const std::string &filename )
  {
#if DUNE_FEM_COMPATIBILITY
    // create read stream 
    XDRReadStream xdr( filename );
    bool success = true;

    // check new version tag
    int newVersionTag = myVersionTag;
    success &= xdr.inout( newVersionTag );
    const bool newVersion = (newVersionTag == myVersionTag);

    // if new version the read type, otherwise newVersionTag is the type info
    int typeVar = (newVersion ? type() : newVersionTag);
    if( newVersion )
      success &= xdr.inout( typeVar );

    // index set type check
    if( (typeVar != 2) && (typeVar != type()) )
    {
      DUNE_THROW( InvalidStateException,
                  "AdaptiveLeafIndexSetBase::read_xdr: wrong type " << typeVar
                  << " given (expected " << type() << ")." );
    }

    if( newVersion )
    {
      // read codim used 
      for( int i = 0; i < numCodimensions; ++i )
        success &= xdr.inout( codimUsed_[ i ] );
    }
    else 
    {
      // it depends on the type whether higher codims were stored
      for( int i = 0; i < numCodimensions; ++i )
        codimUsed_[ i ] = (typeVar == type());
    }

    if( typeVar == type() )
    {
      for( int i = 0; i < numCodimensions; ++i )
      {
        if( codimUsed_[ i ] )
          success &= codimLeafSet_[ i ].processXdr( xdr );
      }
    }
    else
      success &= codimLeafSet_[ 0 ].processXdr( xdr );

    // in parallel runs we have to compress here
    if( grid_.comm().size() > 1 )
      compressed_ = false;
    
    return success;
#else
    // new version using streams 
    try
    {
      XDRFileInStream in( filename );
      return read( in );
    }
    catch( Exception e )
    {
      return false;
    }
#endif
  }


  template< class TraitsImp >
  inline bool AdaptiveLeafIndexSetBase< TraitsImp >
    ::read_xdr ( const std::string &filename, int timestep )
  {
    const char *path = "";
    std::string fnstr = genFilename(path,filename, timestep);
    return read_xdr( fnstr );
  }



  /////////////////////////////////////////////////////////////////////////
  //
  //  --AdaptiveLeafIndexSet 
  //
  /////////////////////////////////////////////////////////////////////////
  template< class GridImp,  
            PartitionIteratorType pt>
  struct AdaptiveLeafIndexSetTraits
  {   
    typedef GridImp GridType;
    static const PartitionIteratorType pitype = pt;
    enum { numCodimensions = GridType :: dimension + 1 };
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
  template < class GridImp,
             PartitionIteratorType pitype = All_Partition > 
  class AdaptiveLeafIndexSet
  : public AdaptiveLeafIndexSetBase< AdaptiveLeafIndexSetTraits< GridImp, pitype > > 
  {
    typedef AdaptiveLeafIndexSetBase< AdaptiveLeafIndexSetTraits< GridImp, pitype > > BaseType;
  public:
    typedef typename BaseType :: GridType GridType;
    //! Constructor
    AdaptiveLeafIndexSet (const GridType & grid) 
      : BaseType(grid) 
    {
    }

    //! return name of index set 
    virtual std::string name () const
    {
      return "AdaptiveLeafIndexSet";
    }
  };

#if ! DUNE_FEM_COMPATIBILITY
  /////////////////////////////////////////////////////////////////////////
  //
  //  --DGAdaptiveLeafIndexSet 
  //
  /////////////////////////////////////////////////////////////////////////
  template< class GridImp,  
            PartitionIteratorType pt>
  struct DGAdaptiveLeafIndexSetTraits
  {   
    typedef GridImp GridType;
    static const PartitionIteratorType pitype = pt;
    // this index set only supports one codimension, codim zero 
    enum { numCodimensions = 1 };
    // first comdimension that is supported (not yet supported)
    enum { startingCodimension = 0 };
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
  template < class GridImp,
             PartitionIteratorType pitype = All_Partition > 
  class DGAdaptiveLeafIndexSet
  : public AdaptiveLeafIndexSetBase< DGAdaptiveLeafIndexSetTraits< GridImp, pitype > > 
  {
    typedef AdaptiveLeafIndexSetBase< DGAdaptiveLeafIndexSetTraits< GridImp, pitype > > BaseType;
  public:
    typedef typename BaseType :: GridType GridType;
    //! Constructor
    DGAdaptiveLeafIndexSet (const GridType & grid) 
      : BaseType(grid) 
    {
    }

    //! return name of index set 
    virtual std::string name () const
    {
      return "DGAdaptiveLeafIndexSet";
    }
  };
#endif


} // end namespace Dune 

#endif
