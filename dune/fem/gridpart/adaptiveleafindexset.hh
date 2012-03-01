#ifndef DUNE_ADAPTIVELEAFINDEXSET_HH
#define DUNE_ADAPTIVELEAFINDEXSET_HH

//- local includes 
#include <dune/common/forloop.hh>

#include <dune/fem/version.hh>
#include <dune/fem/gridpart/dunefemindexsets.hh>
#include <dune/fem/gridpart/codimindexset.hh>

namespace Dune
{
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
  : public ConsecutivePersistentIndexSet< 
            typename TraitsImp :: GridType, 
            AdaptiveIndexSetBase< TraitsImp > 
            >
  {
  protected:  
    typedef typename TraitsImp :: GridPartType GridPartType;
    typedef typename GridPartType :: GridType GridType;
    typedef AdaptiveIndexSetBase< TraitsImp > ThisType;
    typedef ConsecutivePersistentIndexSet< GridType, ThisType > BaseType;

    friend class Conversion< ThisType, Fem :: EmptyIndexSet >;

  public:
    //! dimension of the grid 
    static const int dimension = GridType::dimension;

    //! number of supported codimensions 
    static const int numCodimensions = TraitsImp :: numCodimensions ;

    //! number of supported codimensions 
    static const int intersectionCodimension = TraitsImp :: intersectionCodimension ;

    //! type of index 
    typedef typename BaseType :: IndexType IndexType;

    //! type of codimension 0 Entity (extract from GridPartType)
    typedef typename GridPartType::template Codim< 0 >::EntityType   ElementType;

    //! type of codimension 0 EntityPointer (extract from GridPartType) 
    typedef typename GridPartType::template Codim< 0 >::EntityPointerType  ElementPointerType;

    //! type of intersection iterator 
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType; 

    //! type of intersections
    typedef typename GridPartType :: IntersectionType   IntersectionType; 

  private:
    typedef typename TraitsImp :: CodimIndexSetType  CodimIndexSetType ;

    template< int codim , bool gridHasCodim >
    struct CountElementsBase
    {
      static void apply ( const ThisType &indexSet, const GeometryType &type, int &count )
      {
        if( type.dim() == dimension - codim )
          count = indexSet.template countElements< codim >( type, integral_constant<bool,true>() );
      }
    };

    template< int codim >
    struct CountElementsBase< codim, false >
    {
      static void apply ( const ThisType &indexSet, const GeometryType &type, int &count )
      {
        if( type.dim() == dimension - codim )
          count = indexSet.template countElements< codim >( type, integral_constant<bool,false>() );
      }
    };

    template< int codim >
    struct CountElements  
      : public CountElementsBase< codim, Capabilities :: hasEntity < GridType, codim > :: v >
    {
    };


    template< int codim , bool gridHasCodim >
    struct InsertSubEntitiesBase
    {
      static void apply ( ThisType &indexSet, const ElementType &entity )
      {
        // if codimension is not available return 
        if( ! indexSet.codimAvailable( codim ) ) return ;

        // if codimension is not used return 
        if( !indexSet.codimUsed_[ codim ] ) return;
        
        CodimIndexSetType &codimSet = indexSet.codimLeafSet( codim );

        for( int i = 0; i < entity.template count< codim >(); ++i )
        {
          codimSet.insertSubEntity( entity, i );
        }
      }
    };

    template< int codim >
    struct InsertSubEntitiesBase< codim , false >
    {
      static void apply ( ThisType &indexSet, const ElementType &entity )
      {
      }
    };

    template< int codim >
    struct InsertSubEntities 
      : public InsertSubEntitiesBase< codim , Capabilities :: hasEntity < GridType, codim > :: v >
    {
    };

    template< int codim , bool gridHasCodim >
    struct InsertGhostSubEntitiesBase
    {
      typedef typename GridPartType::template Codim< codim >::EntityPointerType  EntityPointerType;
      typedef typename GridPartType::template Codim< codim >::EntityType         EntityType;

      static void apply ( ThisType &indexSet, const ElementType &entity ,
                          const bool skipGhosts )
      {
        // if codimension is not available return 
        if( ! indexSet.codimAvailable( codim ) ) return ;

        // if codimension is not used return 
        if( !indexSet.codimUsed_[ codim ] ) return;

        CodimIndexSetType &codimSet = indexSet.codimLeafSet( codim );

        for( int i = 0; i < entity.template count< codim >(); ++i )
        {
          EntityPointerType ptr = entity.template subEntity< codim >( i );
          const EntityType &subentity = *ptr;

          if( !skipGhosts || (entity.partitionType() != GhostEntity) )
            codimSet.insertGhost( subentity );
        }
      }
    };

    template< int codim >
    struct InsertGhostSubEntitiesBase< codim, false >
    {
      static void apply ( ThisType &indexSet, const ElementType &entity ,
                          const bool skipGhosts )
      {}
    };

    template< int codim >
    struct InsertGhostSubEntities 
      : public InsertGhostSubEntitiesBase< codim, Capabilities :: hasEntity < GridType, codim > :: v >
    {
    };

    template< int codim , bool gridHasCodim >
    struct CallSetUpCodimSetBase
    {
      static void apply ( const int cd, const ThisType &indexSet )
      {
        // if codimension is not available return 
        if( ! indexSet.codimAvailable( codim ) ) return ;

        if( cd == codim )
          indexSet.template setupCodimSet< codim >(integral_constant<bool,true>());
      }
    };

    template< int codim >
    struct CallSetUpCodimSetBase< codim, false >
    {
      static void apply ( const int cd, const ThisType &indexSet )
      {
        if( cd == codim )
          indexSet.template setupCodimSet< codim >(integral_constant<bool,false>());
      }
    };

    template< int codim >
    struct CallSetUpCodimSet 
      : public CallSetUpCodimSetBase< codim, Capabilities :: hasEntity < GridType, codim > :: v >
    {
    };


    /////////////////////////////////////////////////////
    //  subentity extractor 
    /////////////////////////////////////////////////////

    template < int codim, bool gridHasCodim > 
    struct GetSubEntityBase
    {
      typedef typename GridPartType :: template Codim< codim > :: EntityPointerType  EntityPointer ;
      static EntityPointer subEntity( const ElementType& element, const int subEn ) 
      {
        return element.template subEntity< codim > ( subEn );
      }
    };

    template < int codim > 
    struct GetSubEntityBase< codim, false >
    {
      typedef typename GridPartType :: template Codim< 0 > :: EntityPointerType  EntityPointer ;
      static EntityPointer subEntity( const ElementType& element, const int subEn ) 
      {
        DUNE_THROW(NotImplemented,"stupid grid without entities of codim 1 used");
        return element.template subEntity< 0 > ( 0 );
      }
    };

    struct GetFaceEntity 
      : public GetSubEntityBase< 1, Capabilities :: hasEntity < GridType, 1 > :: v  >
    {
    };

    //! type of codimension 1 entity 
    typedef typename GetFaceEntity :: EntityPointer FacePointerType ;

    //! is true if grid is structured grid 
    enum { StructuredGrid = ! Capabilities::IsUnstructured<GridType>::v };

    // my type, to be revised 
    enum { myType = ( numCodimensions == 1 ) ? ( (StructuredGrid) ? -1 : 665 ) : 6 };

    // max num of codimension (to avoid compiler warnings)
    enum { maxNumCodimension = ((dimension + 1) > numCodimensions) ? dimension + 2 : numCodimensions+1 };

    //! default partition iterator type 
    static const PartitionIteratorType pitype = GridPartType :: indexSetPartitionType ;

    // reference to grid part 
    const GridPartType& gridPart_;
    // Codimension leaf index sets 
    mutable CodimIndexSetType* codimLeafSet_[ numCodimensions ];
    // flag for codim is in use or not 
    mutable bool codimUsed_ [ maxNumCodimension ];

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
    const bool codimAvailable( const int codim ) const 
    {
      return codim < numCodimensions && codim >= 0 ;
    }

    CodimIndexSetType& codimLeafSet( const int codim ) const 
    {   
      assert( codimLeafSet_[ codim ] );
      // assert( codimAvailable( codim ) );
      return *codimLeafSet_[ codim ];
    }

  public:
    //! Constructor
    //AdaptiveIndexSetBase (const GridPartType & gridPart) 
    AdaptiveIndexSetBase (const GridPartType & gridPart)
      : BaseType( gridPart.grid() ) 
      , gridPart_( gridPart )
      , sequence_( dofManager_.sequence() )
      , compressed_(true) // at start the set is compressed 
    {
      // codim 0 is used by default
      codimUsed_[ 0 ] = true;

      // all higher codims are not used by default
      for(int codim = 1; codim < maxNumCodimension; ++codim ) codimUsed_[ codim ] = false ;

      // set the codim of each Codim Set. 
      for(int codim = 0; codim < numCodimensions; ++codim ) 
      {
        if( codim == intersectionCodimension ) 
          codimLeafSet_[ codim ] = new CodimIndexSetType( grid_, 1 );
        else 
          codimLeafSet_[ codim ] = new CodimIndexSetType( grid_, codim );
      }

      // get level-0 view, this is alrady used in GridPtr (DFG parser)
      typedef typename GridType :: LevelGridView MacroViewType;
      MacroViewType macroView = grid_.levelView( 0 );

      // resize vector of geometry types 
      geomTypes_.resize( dimension+1 );
      for(int codim=0; codim <= dimension; ++codim ) 
      {
        // copy geometry types 
        geomTypes_[ codim ] = macroView.indexSet().geomTypes( codim ); 
      }

      // build index set 
      setupIndexSet();
    }

    //! Destructor
    virtual ~AdaptiveIndexSetBase ()
    {
      // set the codim of each Codim Set. 
      for(int codim = 0; codim < numCodimensions; ++codim ) 
      {
        delete codimLeafSet_[ codim ];
        codimLeafSet_[ codim ] = 0;
      }
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
    //! return size of grid entities per level and codim 
    IndexType size ( GeometryType type ) const
    {
      const int codim = dimension - type.dim();

      // use size of codim index set if possible 
      if( codimAvailable( codim ) )
      {
        if( codimUsed_[ codim ] )
          return codimLeafSet( codim ).size();
      }

      assert( geomTypes( codim ).size() == 1 ); 
      int count = 0;
      ForLoop< CountElements, 0, dimension > :: apply( *this, type, count );
      return count;
    }
    
    //! return size of grid entities of given codim 
    IndexType size ( int codim ) const
    {
      assert( codim < numCodimensions );
      if( codim == intersectionCodimension ) 
      {
        return codimLeafSet( codim ).size();
      }
      assert( geomTypes( codim ).size() == 1 ); 
      return size( geomTypes( codim )[0] );
    }
    
    //! returns vector with geometry tpyes this index set has indices for
    const std::vector <GeometryType> & geomTypes (const int codim) const 
    {
      assert( codim >= 0 && codim < int(geomTypes_.size()) );
      return geomTypes_[ codim ];
    }
   
    //! \brief returns true if entity is contained in index set 
    template <class EntityType>
    bool contains (const EntityType & en) const
    {
      enum { codim = EntityType::codimension };
      if( codimAvailable( codim ) )
      {
        assert( codimUsed_[codim] );
        return codimLeafSet( codim ).exists( en ); 
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

    //! make to index numbers consecutive 
    //! return true, if at least one hole existed  
    bool compress ();

  public:
    ////////////////////////////////////////////////////////
    //  index methods 
    //  --index
    ///////////////////////////////////////////////////////
    /* \brief return index for entity */
    template< class Entity >
    IndexType index ( const Entity &entity ) const
    {
      return index< Entity :: codimension >( entity );
    }

    /* \brief return index for entity */
    template< int codim >
    IndexType
    index ( const typename GridType :: template Codim< codim > :: Entity &entity ) const
    {
      if( codimAvailable( codim ) ) 
      {
        if( (codim != 0) && ! codimUsed_[ codim ] )
          setupCodimSet< codim >(integral_constant<bool,true>());

        const CodimIndexSetType &codimSet = codimLeafSet( codim );
        const IndexType idx = codimSet.index( entity );
        assert( (idx >= 0) && (idx < IndexType(codimSet.size())) );
        return idx;
      }
      else 
      {
        DUNE_THROW( NotImplemented, (name() + " does not support indices for codim = ") << codim );
        return -1;
      }
    }

    /* \brief return index for sub entity of given entity and subEntityNumber */
    IndexType
    subIndex ( const ElementType& element,
               int subNumber, unsigned int codim ) const
    {
      if( codimAvailable( codim ) ) 
      {
        if( (codim != 0) && !codimUsed_[ codim ] )
          ForLoop< CallSetUpCodimSet, 0, dimension >::apply( codim, *this );
        
        const CodimIndexSetType &codimSet = codimLeafSet( codim );
        const IndexType idx = codimSet.subIndex( element, subNumber );
        assert( (idx >= 0) && (idx < IndexType(codimSet.size())) );
        return idx;
      }
      else 
      {
        if( (codim != 0) && ! codimUsed_[ codim ] )
        {
          ForLoop< CallSetUpCodimSet, 0, dimension >::apply( codim, *this );
        }
        
        const CodimIndexSetType &codimSet = codimLeafSet( codim );
        const IndexType idx = codimSet.subIndex( element, subNumber );
        assert( (idx >= 0) && (idx < IndexType(codimSet.size())) );
        return idx;
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

        const CodimIndexSetType &codimSet = codimLeafSet( codim );

        // get corresponding face entity pointer 
        FacePointerType face = getIntersectionFace( intersection );
        
        const IndexType idx = codimSet.index( *face );
        assert( (idx >= 0) && (idx < IndexType(codimSet.size())) );
        return idx;
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
    //////////////////////////////////////////////////////////
    //
    //  DoF adjustment methods
    //
    //////////////////////////////////////////////////////////

    //! return number of holes of the sets indices 
    int numberOfHoles ( const int codim ) const
    {
      if( codimAvailable( codim ) ) 
      {
        assert( codimUsed_[codim] );
        return codimLeafSet( codim ).numberOfHoles(); 
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
        return codimLeafSet( codim ).oldIndex( hole ); 
      }
      else 
      {
        DUNE_THROW( NotImplemented, (name() + " does not support indices for codim = ") << codim );
        return -1; 
      }
    }

    //! return new index, for dof manager only returns index 
    int newIndex (const int hole , const int codim ) const
    {
      if( codimAvailable( codim ) ) 
      {
        assert( codimUsed_[codim] );
        return codimLeafSet( codim ).newIndex( hole ); 
      }
      else 
      {
        DUNE_THROW( NotImplemented, (name() + " does not support indices for codim = ") << codim );
        return -1; 
      }
    }

  protected:
    // memorise index 
    void insertIndex ( const ElementType &entity );

    // memorise indices for all intersections  
    void insertIntersections ( const ElementType &entity ) const;

    // insert index temporarily
    void insertTemporary ( const ElementType &entity );

    // set indices to unsed so that they are cleaned on compress
    void removeIndex ( const ElementType &entity );

    // check whether entity can be inserted or not 
    void checkHierarchy ( const ElementType &entity, const bool wasNew );

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
    void setupCodimSet (const integral_constant<bool,true> &hasEntities) const;
    template< int codim >
    void setupCodimSet (const integral_constant<bool,false> &hasEntities) const;

    // mark indices that are still used (and give new indices to new intersections)
    void setupIntersections () const;

    // count elements by iterating over grid and compare 
    // entities of given codim with given type 
    template< int codim >
    inline int countElements ( GeometryType type, const integral_constant<bool,true> &hasEntities ) const;
    template< int codim >
    inline int countElements ( GeometryType type, const integral_constant<bool,false> &hasEntities ) const;
    
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

  protected:
    FacePointerType getIntersectionFace( const IntersectionType& intersection ) const 
    {
      ElementPointerType inside = intersection.inside();
      return getIntersectionFace( intersection, *inside );
    }

    FacePointerType getIntersectionFace( const IntersectionType& intersection, 
                                         const ElementType& inside ) const 
    {
      if( ! intersection.conforming() && intersection.neighbor() ) 
      {
        const ElementPointerType outsideEp = intersection.outside();
        const ElementType& outside = *outsideEp ;
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
        if( codimUsed_[ codim ] )
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
      codimLeafSet( codim ).clearHoles();

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
      codimLeafSet( codim ).prepareCompress();
    }

    // mark all indices still needed 
    setupIndexSet();

    // true if a least one index is moved
    bool haveToCopy = codimLeafSet( 0 ).compress();
    for( int codim = 1; codim < numCodimensions; ++codim )
    {
      if( codimUsed_[ codim ] ) 
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
  AdaptiveIndexSetBase< TraitsImp >::insertIndex ( const ElementType &entity )
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
        ForLoop< InsertGhostSubEntities, 1, dimension >::apply( *this, entity, skipGhosts );
    }
    else 
#endif // HAVE_MPI
    {
      codimLeafSet( 0 ).insert( entity );
      // only for index sets upporting more than one codim 
      if( numCodimensions > 1 )
        ForLoop< InsertSubEntities, 1, dimension >::apply( *this, entity );

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
  AdaptiveIndexSetBase< TraitsImp >::insertIntersections ( const ElementType &entity ) const
  {
    codimLeafSet( intersectionCodimension ).resize();

    const IntersectionIteratorType endiit = gridPart_.iend( entity );
    for( IntersectionIteratorType iit = gridPart_.ibegin( entity ); 
         iit != endiit ; ++ iit ) 
    {
      // get intersection 
      const IntersectionType& intersection = *iit ;

      // get correct face pointer 
      FacePointerType face = getIntersectionFace( intersection, entity );

      // insert face into index set 
      codimLeafSet( intersectionCodimension ).insert( *face );
    }
  }

  template< class TraitsImp >
  inline void
  AdaptiveIndexSetBase< TraitsImp >::insertTemporary( const ElementType &entity )
  {
    insertIndex( entity );
    codimLeafSet( 0 ).markForRemoval( entity );
  }

  template< class TraitsImp >
  inline void
  AdaptiveIndexSetBase< TraitsImp >::removeIndex( const ElementType &entity )
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
      if( codimUsed_[ codim ] )
        codimLeafSet( codim ).resetUsed();
    }

    typedef typename GridPartType 
      ::template Codim< 0 > :: template Partition< pt > :: IteratorType  Iterator;

    const Iterator end  = gridPart_.template end< 0, pt >();
    for( Iterator it = gridPart_.template begin< 0, pt >(); it != end; ++it )
      insertIndex( *it );
  }

  template< class TraitsImp >
  inline void
  AdaptiveIndexSetBase< TraitsImp >::clear()
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
    // so we can use SGrid for the visualization :(
    if( StructuredGrid  &&
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
      if( codimUsed_[ codim ] )
        codimLeafSet( codim ).resetUsed(); 
    }
    
    // get macro iterator 
    typedef typename GridType
      ::template Codim< 0 >::template Partition< pt >::LevelIterator
      Iterator;

    typedef typename GridType :: template Partition< All_Partition > :: LevelGridView LevelGridView ;
    LevelGridView macroView = grid_.levelView( 0 );

    const Iterator macroend = macroView.template end< 0, pt >();
    for( Iterator macroit = macroView.template begin< 0, pt >(); 
         macroit != macroend; ++macroit )
      checkHierarchy( *macroit, false );
  }


  template< class TraitsImp >
  template< int codim >
  inline void
  AdaptiveIndexSetBase< TraitsImp >::setupCodimSet (const integral_constant<bool,true>&) const
  {
    // if codim is not available do nothing 
    if( ! codimAvailable( codim ) ) return ;

    // resize if necessary 
    codimLeafSet( codim ).resize();
    
    // walk over grid parts entity set and insert entities
    typedef typename GridPartType
      ::template Codim< codim >::template Partition< pitype > :: IteratorType Iterator;

    const Iterator end = gridPart_.template end< codim, pitype >();
    for( Iterator it = gridPart_.template begin< codim, pitype >(); it != end; ++it )
      codimLeafSet( codim ).insert( *it );

    // mark codimension as used
    codimUsed_[ codim ] = true;
  }
  template< class TraitsImp >
  template< int codim >
  inline void
  AdaptiveIndexSetBase< TraitsImp >::setupCodimSet (const integral_constant<bool,false>&) const
  {
    // if codim is not available do nothing 
    if( ! codimAvailable( codim ) ) return ;

    // resize if necessary 
    codimLeafSet( codim ).resize();

    typedef typename GridPartType
      ::template Codim< 0 >::template Partition< pitype > :: IteratorType Iterator;

    const Iterator end = gridPart_.template end< 0, pitype >();
    for( Iterator it = gridPart_.template begin< 0, pitype >(); it != end; ++it )
    {
      for (int i=0;i<it->template count<codim>();++i)
      {
        if (! codimLeafSet( codim ).exists( *it, i) )
          codimLeafSet( codim ).insertSubEntity( *it,i );
      }
    }

    // mark codimension as used
    codimUsed_[ codim ] = true;
  }


  template< class TraitsImp >
  inline void
  AdaptiveIndexSetBase< TraitsImp >::setupIntersections() const
  {
    // if intersectionCodimension < 0 then this feature is disabled 
    if( intersectionCodimension < 0 ) return ;

    // do nothing if insections are already available 
    if( codimUsed_[ intersectionCodimension ] ) return ;

    // resize if necessary 
    codimLeafSet( intersectionCodimension ).resize();
    
    // walk over grid parts entity set and insert entities
    typedef typename GridPartType
      ::template Codim< 0 >::template Partition< pitype > :: IteratorType Iterator;

    const Iterator end = gridPart_.template end< 0, pitype >();
    for( Iterator it = gridPart_.template begin< 0, pitype >(); it != end; ++it )
    {
      // insert all intersections of this entity 
      insertIntersections( *it );
    }

    // mark codimension as used
    codimUsed_[ intersectionCodimension ] = true;
  }

  template< class TraitsImp >
  template< int codim >
  inline int
  AdaptiveIndexSetBase< TraitsImp >::countElements ( GeometryType type, const integral_constant<bool,true>& ) const
  {
    typedef typename GridPartType
      ::template Codim< codim > :: template Partition< pitype > :: IteratorType Iterator;

    const Iterator begin = gridPart_.template begin< codim, pitype >();
    const Iterator end = gridPart_.template end< codim, pitype >();
    int count = 0;
    for( Iterator it = begin; it != end; ++it )
    {
      if( it->type() == type )
        ++count;
    }
    return count; 
  }
  template< class TraitsImp >
  template< int codim >
  inline int
  AdaptiveIndexSetBase< TraitsImp >::countElements ( GeometryType type, const integral_constant<bool,false>& ) const
  {
    // make sure codimension is enabled 
    assert( codimAvailable( codim ) );

    // resize if necessary 
    codimLeafSet( codim ).resize();

    typedef typename GridPartType
      ::template Codim< 0 >::template Partition< pitype > :: IteratorType Iterator;

    const Iterator end = gridPart_.template end< 0, pitype >();
    int count = 0;
    for( Iterator it = gridPart_.template begin< 0, pitype >(); it != end; ++it )
    {
      for (int i=0;i<it->template count<codim>();++i)
      {
        if (! codimLeafSet( codim ).exists( *it, i) )
        {
          codimLeafSet( codim ).insertSubEntity( *it,i );
          if ( Dune::GenericReferenceElements< typename GridPartType::ctype, GridPartType::dimension >::
             general( it->type() ).type( i,codim ) == type )
            ++count;
        }
      }
    }

    // mark codimension as used
    codimUsed_[ codim ] = true;

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
      out << codimUsed_[ i ];

    // write all sets 
    for( int i = 0; i < numCodimensions; ++i )
    {
      if( codimUsed_[ i ] )
        codimLeafSet( i ).write( out );
    }
    
    // if we got until here writing was sucessful
    return true;
  }

  template< class TraitsImp >
  inline bool AdaptiveIndexSetBase< TraitsImp >
    ::write_xdr ( const std::string &filename ) const
  {
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
  }


  template< class TraitsImp >
  inline bool AdaptiveIndexSetBase< TraitsImp >
    ::write_xdr ( const std::string &filename, int timestep ) const
  {
    const char *path = "";
    std::string fnstr = genFilename( path, filename, timestep );
    return write_xdr( fnstr );
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
        DUNE_THROW( InvalidStateException,
                    "AdaptiveIndexSetBase::read: got " << storedName 
                    << " (expected " << storedName << ")." );
      }
    }

    // read number of codimensions 
    int numCodim;
    in >> numCodim; 

    // make sure everything is correct 
    assert( numCodim == numCodimensions );

    // read codim used 
    for( int i = 0; i < numCodimensions; ++i )
      in >> codimUsed_[ i ];

    for( int i = 0; i < numCodimensions; ++i )
    {
      if( codimUsed_[ i ] )
        codimLeafSet( i ).read( in );
    }

    // in parallel runs we have to compress here
    if( grid_.comm().size() > 1 )
      compressed_ = false;
    
    // if we got until here reading was sucessful 
    return true;
  }

  template< class TraitsImp >
  inline bool AdaptiveIndexSetBase< TraitsImp >
    ::read_xdr ( const std::string &filename )
  {
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
  }


  template< class TraitsImp >
  inline bool AdaptiveIndexSetBase< TraitsImp >
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
  template< class GridPartImp >
  struct AdaptiveLeafIndexSetTraits
  {   
    // type of grid part  
    typedef GridPartImp GridPartType;
    // type of grid 
    typedef typename GridPartType :: GridType GridType;
    // number of codimensions 
    enum { numCodimensions = GridType :: dimension + 1 };
    // first comdimension that is supported (not yet supported)
    enum { startingCodimension = 0 };
    // intersection codimensions (this is usually dimension + 1 )
    enum { intersectionCodimension = -1 };
    // type of codimension index set  
    typedef CodimIndexSet< GridType >  CodimIndexSetType; 
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
  public:
    typedef typename BaseType :: GridPartType GridPartType;
    //! Constructor
    AdaptiveLeafIndexSet (const GridPartType & gridPart) 
      : BaseType(gridPart) 
    {
    }

    //! return name of index set 
    virtual std::string name () const
    {
      return "AdaptiveLeafIndexSet";
    }
  };


  /////////////////////////////////////////////////////////////////////////
  //
  //  --IntersectionAdaptiveLeafIndexSet 
  //
  /////////////////////////////////////////////////////////////////////////
  template< class GridPartImp >
  struct IntersectionAdaptiveLeafIndexSetTraits
  {   
    // type of grid part  
    typedef GridPartImp GridPartType;
    // type of grid 
    typedef typename GridPartType :: GridType GridType;
    // number of codimensions 
    enum { numCodimensions = GridType :: dimension + 2 };
    // intersection codimensions (this is usually dimension + 1 )
    enum { intersectionCodimension = numCodimensions - 1 };
    // first comdimension that is supported (not yet supported)
    enum { startingCodimension = 0 };
    // type of codimension index set  
    typedef CodimIndexSet< GridType >  CodimIndexSetType; 
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
  public:
    typedef typename BaseType :: GridPartType GridPartType;
    //! Constructor
    IntersectionAdaptiveLeafIndexSet (const GridPartType & gridPart) 
      : BaseType(gridPart) 
    {
    }

    //! return name of index set 
    virtual std::string name () const
    {
      return "IntersectionAdaptiveLeafIndexSet";
    }
  };

  /////////////////////////////////////////////////////////////////////////
  //
  //  --DGAdaptiveLeafIndexSet 
  //
  /////////////////////////////////////////////////////////////////////////
  template< class GridPartImp >
  struct DGAdaptiveLeafIndexSetTraits
  {   
    // type of grid part 
    typedef GridPartImp GridPartType;
    // type of grid 
    typedef typename GridPartType :: GridType GridType;
    // this index set only supports one codimension, codim zero 
    enum { numCodimensions = 1 };
    // first comdimension that is supported (not yet supported)
    enum { startingCodimension = 0 };
    // intersection codimensions (this is usually dimension + 1 )
    enum { intersectionCodimension = -1 };
    // type of codimension index set  
    typedef CodimIndexSet< GridType >  CodimIndexSetType; 
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
  public:
    typedef typename BaseType :: GridPartType GridPartType;
    //! Constructor
    DGAdaptiveLeafIndexSet (const GridPartType & gridPart) 
      : BaseType(gridPart) 
    {
    }

    //! return name of index set 
    virtual std::string name () const
    {
      return "DGAdaptiveLeafIndexSet";
    }
  };

  /////////////////////////////////////////////////////////////////////////
  //
  //  --IdBasedLeafIndexSet 
  //
  /////////////////////////////////////////////////////////////////////////
  template< class GridPartImp >
  struct IdBasedLeafIndexSetTraits
  {   
    // type of grid part 
    typedef GridPartImp GridPartType;
    // type of grid 
    typedef typename GridPartType :: GridType GridType;
    // this index set only supports one codimension, codim zero 
    enum { numCodimensions = GridType :: dimension + 1 };
    // first comdimension that is supported (not yet supported)
    enum { startingCodimension = 0 };
    // intersection codimensions (this is usually dimension + 1 )
    enum { intersectionCodimension = -1 };
    // type of codimension index set  
    typedef CodimIndexSet< GridType >  CodimIndexSetType; 
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
  class IdBasedLeafIndexSet
  : public AdaptiveIndexSetBase< IdBasedLeafIndexSetTraits< GridPartImp > > 
  {
    typedef AdaptiveIndexSetBase< IdBasedLeafIndexSetTraits< GridPartImp > > BaseType;
  public:
    typedef typename BaseType :: GridPartType GridPartType;
    //! Constructor
    DUNE_VERSION_DEPRECATED(1,3,remove) 
    IdBasedLeafIndexSet (const GridPartType & gridPart) 
      : BaseType(gridPart) 
    {
    }

    //! return name of index set 
    virtual std::string name () const
    {
      return "IdBasedLeafIndexSet";
    }
  };

} // end namespace Dune 

#endif
