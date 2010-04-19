#ifndef DUNE_FILTEREDINDEXSET_HH
#define DUNE_FILTEREDINDEXSET_HH

#warning "Deprecated header, don't use this index set anymore!"

//- local includes 
#include <dune/common/forloop.hh>
#include <dune/fem/gridpart/dunefemindexsets.hh>
#include <dune/fem/gridpart/codimindexset.hh>
#include <dune/fem/io/file/xdrio.hh>

namespace Dune
{
  template< class GridPartImp, class FilterImp, bool >
  class FilteredGridPart;

  /** \class FilteredIndexSet
   *  \brief consecutive, persistent index set for the grid part given 
   *
   *  This index set generates a consecutive index out of the unique global
   *  index of each entity. It can be used instead of the default grid index sets
   *  and can be generated for each grid implementation.
   *
   */
  template <class GridPartImp, class FilterImp >
  class FilteredIndexSet
  : public ConsecutiveIndexSet< typename GridPartImp :: GridType, 
                                FilteredIndexSet< GridPartImp, FilterImp > >
  {
    typedef typename GridPartImp :: GridType GridType;
    typedef FilteredGridPart< GridPartImp, FilterImp , true > GridPartType;

    typedef FilteredIndexSet< GridPartImp, FilterImp > ThisType;
    typedef ConsecutiveIndexSet< GridType, ThisType > BaseType;

    friend class Conversion< ThisType, EmptyIndexSet >;

  public:
    static const int dimension = GridType::dimension;

    static const int ncodim = dimension + 1;

    static const PartitionIteratorType pitype = GridPartImp :: indexSetPartitionType;

    //! type of index 
    typedef typename BaseType :: IndexType IndexType;

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
      static void apply ( const ThisType &indexSet, const ElementType &entity )
      {
        if( ! indexSet.codimUsed_[ codim ] )
          return;
        
        const HIndexSetType &hIndexSet = indexSet.hIndexSet_;
        CodimIndexSetType &codimSet = indexSet.codimLeafSet_[ codim ];

        const int subEntities = entity.template count< codim >();
        for( int i = 0; i < subEntities ; ++i )
        {
          codimSet.insert( hIndexSet.subIndex( entity, i, codim ) );
        }
      }
    };

    template< int codim >
    struct CallSetUpCodimSet
    {
      static void apply ( const int cd, const ThisType &indexSet )
      {
        if( cd == codim )
          indexSet.template setUpCodimSet< codim >();
      }
    };

    // my type, to be revised 
    enum { myType = 6 };
    enum { myVersionTag = -665 };

    const GridPartType& gridPart_;
    const HIndexSetType &hIndexSet_;
    mutable CodimIndexSetType codimLeafSet_[ ncodim ];
    // flag for codim is in use or not 
    mutable bool codimUsed_ [ncodim];
    
    //! flag is tru if set is in compressed status
    mutable bool compressed_;

  protected:
    using BaseType::grid_;

  public:
    //! type traits of this class (see defaultindexsets.hh)
    typedef DefaultLeafIteratorTypes<GridType> Traits; 

    //! Constructor
    DUNE_VERSION_DEPRECATED(1,2,remove) 
    FilteredIndexSet (const GridPartType & gridPart) 
      : BaseType(gridPart.grid()) 
      , gridPart_( gridPart )  
      , hIndexSet_( SelectorType::hierarchicIndexSet(gridPart.grid()) ) 
      , compressed_(true) // at start the set is compressed 
    {
      // codim 0 is used by default
      codimUsed_[0] = true;

      // all higher codims are not used by default
      for(int i=1; i<ncodim; ++i) codimUsed_[i] = false ;
      
      // set the codim of each Codim Set. 
      for(int i=0; i<ncodim; ++i) codimLeafSet_[i].setCodim( i );

      resizeVectors();
      // give all entities that lie below the old entities new numbers 
      markAllUsed ();
    }

    //! Destructor
    virtual ~FilteredIndexSet ()
    {}

    //! return type of index set, for GrapeDataIO
    int type () const
    {
      return myType;
    }

    //! return name of index set 
    std::string name () const
    {
      return "FilteredIndexSet";
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
      if( codimUsed_[ codim ] )
        return codimLeafSet_[ codim ].size();

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
      assert( codimUsed_[codim] );
      return codimLeafSet_[codim].exists( hIndexSet_.index( en ) ); 
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
      markAllUsed();
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
      if( (codim != 0) && !codimUsed_[ codim ] )
        setUpCodimSet< codim >();

      //assertCodimSetSize( codim );
      const CodimIndexSetType &codimSet = codimLeafSet_[ codim ];
      const int hIdx = hIndexSet_.template index( entity );
      const int idx = codimSet.index( hIdx );
      assert( (idx >= 0) && (idx < codimSet.size()) );
      return idx;
    }

    template< int codim >
    IndexType DUNE_DEPRECATED
    subIndex ( const typename GridType :: template Codim< 0 > :: Entity &entity,
               int subNumber ) const
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

    IndexType
    subIndex ( const typename GridType::template Codim< 0 >::Entity &entity,
               int subNumber, unsigned int codim ) const
    {
      if( (codim != 0) && !codimUsed_[ codim ] )
        ForLoop< CallSetUpCodimSet, 0, dimension >::apply( codim, *this );
      
      const CodimIndexSetType &codimSet = codimLeafSet_[ codim ];
      const int hIdx = hIndexSet_.subIndex( entity, subNumber, codim );
      const int idx = codimSet.index( hIdx );
      assert( (idx >= 0) && (idx < codimSet.size()) );
      return idx;
    }

    //! return number of holes of the sets indices 
    int numberOfHoles ( const int codim ) const
    {
      assert( codimUsed_[codim] );
      return codimLeafSet_[codim].numberOfHoles(); 
    }

    //! return old index, for dof manager only 
    int oldIndex (const int hole, const int codim ) const
    {
      assert( codimUsed_[codim] );
      return codimLeafSet_[codim].oldIndex( hole ); 
    }

    //! return new index, for dof manager only returns index 
    int newIndex (const int hole , const int codim ) const
    {
      assert( codimUsed_[codim] );
      return codimLeafSet_[codim].newIndex( hole ); 
    }

  protected:
    // memorise index 
    void insertIndex ( const ElementType &entity );

    // set indices to unsed so that they are cleaned on compress
    void removeIndex ( const ElementType &entity );

    // mark indices that are still used (and give new indices to new elements)
    void markAllUsed (); 

    // mark indices that are still used (and give new indices to new elements)
    template< int codim >
    void setUpCodimSet () const;

    // count elements by iterating over grid and compare 
    // entities of given codim with given type 
    template< int codim >
    inline int countElements ( GeometryType type ) const;
    
  public:
    //! write indexset to xdr file 
    bool write_xdr( const std::string &filename );
    //! write indexset to xdr file 
    bool write_xdr( const std::string &filename, int timestep );

    //! read index set from given xdr file 
    bool read_xdr( const std::string &filename, int timestep );
    //! read index set from given xdr file 
    bool read_xdr( const std::string &filename );
  };



  template< class GridPartType, class FilterImp >
  inline void
  FilteredIndexSet< GridPartType, FilterImp >::resizeVectors ()
  {
    codimLeafSet_[ 0 ].resize( hIndexSet_.size( 0 ) );

    for( int codim = 1; codim < ncodim; ++codim )
    {
      if( codimUsed_[ codim ] )
        codimLeafSet_[ codim ].resize( hIndexSet_.size( codim ) );
    }
  }


  template< class GridPartType, class FilterImp >
  inline bool
  FilteredIndexSet< GridPartType, FilterImp >::compress ()
  {
    // reset list of holes in any case
    for( int i = 0; i < ncodim; ++i )
      codimLeafSet_[ i ].clearHoles();

    if( compressed_ )
    {
      return false;
    }

    // iterate over the leaf level and mark all leaf entities used
    markAllUsed();

    // true if a least one index is moved
    bool haveToCopy = codimLeafSet_[ 0 ].compress();
    for( int i = 1; i < ncodim; ++i )
    {
      if( codimUsed_[i] ) 
        haveToCopy |= codimLeafSet_[ i ].compress();
    }

    // now status is compressed
    compressed_ = true;

    return haveToCopy;
  }


  template< class GridPartType, class FilterImp >
  inline void
  FilteredIndexSet< GridPartType, FilterImp >::insertIndex ( const ElementType &entity )
  {
    const int index = hIndexSet_.index( entity );
    codimLeafSet_[ 0 ].insert( index );
    ForLoop< InsertSubEntities, 1, dimension >::apply( *this, entity );

    assert( codimLeafSet_[ 0 ].exists( index ) );

    // now consecutivity is no longer guaranteed
    compressed_ = false;
  }


  template< class GridPartType, class FilterImp >
  inline void
  FilteredIndexSet< GridPartType, FilterImp >::removeIndex( const ElementType &entity )
  {
    // remove entities (only mark them as unused)
    codimLeafSet_[ 0 ].remove( hIndexSet_.index( entity ) );

    // don't remove higher codim indices (will be done on compression

    // now consecutivity is no longer guaranteed
    compressed_ = false;
  }

  template< class GridPartType, class FilterImp >
  inline void
  FilteredIndexSet< GridPartType, FilterImp >::markAllUsed ()
  {
    // make correct size of vectors
    resizeVectors();

    // mark all indices as unused
    for( int i = 0; i < ncodim; ++i )
    {
      if( codimUsed_[ i ] )
        codimLeafSet_[ i ].set2Unused();
    }
    
    // walk over leaf level on locate all needed entities  
    const PartitionIteratorType pt = All_Partition;
    //const PartitionIteratorType pt = pitype;

    typedef typename GridPartType
      ::template Codim< 0 >::template Partition< pt >:: IteratorType Iterator;

    const Iterator end  = gridPart_.template end< 0, pt >();
    for( Iterator it = gridPart_.template begin< 0, pt >(); it != end; ++it )
      insertIndex( *it );
  }

  template< class GridPartType, class FilterImp >
  template< int codim >
  inline void
  FilteredIndexSet< GridPartType, FilterImp >::setUpCodimSet () const
  {
    // resize if necessary 
    codimLeafSet_[ codim ].resize( hIndexSet_.size( codim ) );
    
    // mark codimension as used
    codimUsed_[codim] = true;

    // walk over leaf level on locate all needed entities  
    typedef typename GridPartType :: template Codim< 0 >:: 
      template Partition< pitype >:: IteratorType Iterator;

    const Iterator end = gridPart_.template end< 0, pitype >();
    for( Iterator it = gridPart_.template begin< 0, pitype >(); 
      it != end; ++it )
    {
      const ElementType& entity = *it;

      if( codim == 0 ) 
        codimLeafSet_[ codim ].insert( hIndexSet_.index ( entity ));
      else 
      {
        InsertSubEntities< codim > :: apply( *this, entity );
      }
    }
  }


  template< class GridPartType, class FilterImp >
  template< int cd >
  inline int
  FilteredIndexSet< GridPartType, FilterImp >::countElements ( GeometryType type ) const
  {
    if( cd == 0 )
    {
      enum { codim = 0 };
      typedef typename GridPartType
        ::template Codim< codim >::template Partition< pitype >:: IteratorType
        Iterator;

      int count = 0;
      const Iterator end = gridPart_.template end< codim, pitype >();
      for( Iterator it = gridPart_.template begin< codim, pitype >();
           it != end; ++it )
      {
        if( it->type() == type )
          ++count;
      }
      return count; 
    }
    else 
    {
      if( !codimUsed_[ cd ] )
        setUpCodimSet< cd >();

      return codimLeafSet_[ cd ].size();
    }
  }


  template< class GridPartType, class FilterImp >
  inline bool FilteredIndexSet< GridPartType, FilterImp >::
  write_xdr ( const std::string &filename )
  {
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
    for( int i = 0; i < ncodim; ++i )
      success &= xdr.inout( codimUsed_[ i ] );

    // write all sets 
    for( int i = 0; i < ncodim; ++i )
      success &= codimLeafSet_[ i ].processXdr( xdr );
    
    return success;
  }


  template< class GridPartType, class FilterImp >
  inline bool FilteredIndexSet< GridPartType, FilterImp >::
  write_xdr ( const std::string &filename, int timestep )
  {
    const char *path = "";
    std::string fnstr = genFilename( path, filename, timestep );
    return write_xdr( fnstr );
  }


  template< class GridPartType, class FilterImp >
  inline bool FilteredIndexSet< GridPartType, FilterImp >::
  read_xdr ( const std::string &filename )
  {
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
                  "FilteredIndexSet::read_xdr: wrong type " << typeVar
                  << " given (expected " << type() << ")." );
    }

    if( newVersion )
    {
      // read codim used 
      for( int i = 0; i < ncodim; ++i )
        success &= xdr.inout( codimUsed_[ i ] );
    }
    else 
    {
      // it depends on the type whether higher codims were stored
      for( int i = 0; i < ncodim; ++i )
        codimUsed_[ i ] = (typeVar == type());
    }

    if( typeVar == type() )
    {
      for( int i = 0; i < ncodim; ++i )
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
  }


  template< class GridPartType, class FilterImp >
  inline bool FilteredIndexSet< GridPartType, FilterImp >::
  read_xdr ( const std::string &filename, int timestep )
  {
    const char *path = "";
    std::string fnstr = genFilename(path,filename, timestep);
    return read_xdr( fnstr );
  }

} // end namespace Dune 

#endif
