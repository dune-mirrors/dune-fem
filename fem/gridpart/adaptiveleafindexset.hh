#ifndef DUNE_ADAPTIVELEAFINDEXSET_HH
#define DUNE_ADAPTIVELEAFINDEXSET_HH

//- local includes 
#include <dune/fem/misc/forloop.hh>
#include <dune/fem/gridpart/dunefemindexsets.hh>
#include <dune/fem/gridpart/codimindexset.hh>
#include <dune/fem/io/file/xdrio.hh>

namespace Dune
{

/** \class AdaptiveLeafIndexSet
 *  \brief consecutive, persistent index set for the leaf level based on the
 *         grid's hierarchy index set
 *
 *  This index set generates a consecutive leaf index out of the unique global
 *  index of each entity. It can be used instead of the default grid index sets
 *  and can be generated for each grid implementation.
 *
 *  \note Only codimension 0 is working at the moment.
 *
 *  \todo Support higher codimensions
 */
template <class GridType, PartitionIteratorType pitype = All_Partition >
class AdaptiveLeafIndexSet
: public ConsecutivePersistentIndexSet
  < GridType, AdaptiveLeafIndexSet< GridType, pitype >, DefaultLeafIteratorTypes< GridType > >
{
  typedef AdaptiveLeafIndexSet< GridType, pitype > ThisType;
  typedef ConsecutivePersistentIndexSet
    < GridType, ThisType, DefaultLeafIteratorTypes< GridType > >
    BaseType;

  friend class Conversion< ThisType, EmptyIndexSet >;

public:
  static const int dimension = GridType :: dimension;

  static const int ncodim = dimension + 1;

  //! type of index 
  typedef typename BaseType :: IndexType IndexType;

  typedef typename GridType :: template Codim< 0 > :: Entity ElementType;

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
      if( !indexSet.codimUsed_[ codim ] )
        return;
      
      const HIndexSetType &hIndexSet = indexSet.hIndexSet_;
      CodimIndexSetType &codimSet = indexSet.codimLeafSet_[ codim ];
      for( int i = 0; i < entity.template count< codim >(); ++i )
        codimSet.insert( hIndexSet.template subIndex< codim >( entity, i ) );
    }
  };

  template< int codim >
  struct InsertGhostSubEntities
  {
    static void apply ( ThisType &indexSet, const ElementType &entity ,
                        const bool skipGhosts )
    {
      if( !indexSet.codimUsed_[ codim ] )
        return;
      
      typedef typename GridType :: template Codim<codim> :: EntityPointer EntityPointerType;
      typedef typename GridType :: template Codim<codim> :: Entity EntityCodimType;

      const HIndexSetType &hIndexSet = indexSet.hIndexSet_;
      CodimIndexSetType &codimSet = indexSet.codimLeafSet_[ codim ];

      for( int i = 0; i < entity.template count< codim >(); ++i )
      {
        // get entity point to check partition type 
        EntityPointerType ep = entity.template entity< codim > (i);
        const EntityCodimType& e = *ep;
        if( e.partitionType() == GhostEntity )
        {
          if( ! skipGhosts )
          {
            // insert ghost entity 
            codimSet.insertGhost( hIndexSet.index( e ) );
          }
        }
        else
        {
          // insert border entity  
          codimSet.insertGhost( hIndexSet.index( e ) );
        }
      }
    }
  };

  // my type, to be revised 
  enum { myType = 6 };
  enum { myVersionTag = -665 };

  const HIndexSetType &hIndexSet_;
  mutable CodimIndexSetType codimLeafSet_[ ncodim ];
  // flag for codim is in use or not 
  mutable bool codimUsed_ [ncodim];
  
  // true if any of the higher codims is used 
  mutable bool higherCodims_;

  //! flag is tru if set is in compressed status
  mutable bool compressed_;

  // actual sequence number 
  int sequence_;

public:
  //! type traits of this class (see defaultindexsets.hh)
  typedef DefaultLeafIteratorTypes<GridType> Traits; 

  //! Constructor
  AdaptiveLeafIndexSet (const GridType & grid) 
    : BaseType(grid) 
    , hIndexSet_( SelectorType::hierarchicIndexSet(grid) ) 
    , higherCodims_( false ) // higherCodims are not used by default 
    , compressed_(true) // at start the set is compressed 
    , sequence_( this->dofManager_.sequence() )
  {
    // codim 0 is used by default
    codimUsed_[0] = true;

    // all higher codims are not used by default
    for(int i=1; i<ncodim; ++i) codimUsed_[i] = higherCodims_;
    
    // set the codim of each Codim Set. 
    for(int i=0; i<ncodim; ++i) codimLeafSet_[i].setCodim( i );

    resizeVectors();
    // give all entities that lie below the old entities new numbers 
    markAllUsed ();
  }

  //! Destructor
  virtual ~AdaptiveLeafIndexSet () {};

  //! return type of index set, for GrapeDataIO
  int type () const { return myType; }

  //! return name of index set 
  std::string name () const { return "AdaptiveLeafIndexSet"; }

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
 
#if defined INDEXSET_HAS_ITERATORS
  /** @brief Iterator to one past the last entity of given codim for partition type
   *  Here the grids leaf iterator is used 
   */
  template<int cd, PartitionIteratorType pt>
  typename DefaultLeafIteratorTypes<GridType>::template Codim<cd>::
    template Partition<pt>::Iterator end () const
  {
    return this->grid_.template leafend<cd,pt> ();
  }

  /** @brief Iterator to first entity of given codimension and partition type.
   *  Here the grids leaf iterator is used 
   */
  template<int cd, PartitionIteratorType pt>
  typename DefaultLeafIteratorTypes<GridType>::template Codim<cd>::
    template Partition<pt>::Iterator begin () const
  {
    return this->grid_.template leafbegin<cd,pt> ();
  }
#endif
 
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
  void insertEntity(const typename GridType::template Codim<0>::Entity & en )  
  {
    // here we have to add the support of higher codims 
    resizeVectors();
    
    this->insertIndex( en );
  }

  //! Unregister entity from set 
  void removeEntity(const typename GridType::template Codim<0>::Entity & en )
  {
    this->removeIndex( en ); 
  }

  //! reallocate the vector for new size
  void resizeVectors()
  {
    codimLeafSet_[0].resize( hIndexSet_.size(0) );
    if(higherCodims_)
    {
      for(int i=1; i<ncodim; ++i) 
      {  
        if(codimUsed_[i]) 
        {
          codimLeafSet_[i].resize( hIndexSet_.size(i) );
        }
      }
    }
  }

  //! if grid has changed, resize index vectors, and create 
  //! indices for new entities, new entities are entities that 
  //! lie below the old entities 
  //- --resize 
  void resize () 
  {
    //std::cout << "Resizing the index set " << this << " \n"; 
    resizeVectors();

    // give all entities that lie below the old entities new numbers 
    markAllBelowOld ();
  }

  //! make to index numbers consecutive 
  //! return true, if at least one hole existed  
  // --compress 
  bool compress ()
  {
    for(int i=0; i<ncodim; ++i) 
    {
      // reset list of holes in any case 
      codimLeafSet_[i].clearHoles();
    }

    // if set already compress, do noting (only for serial runs) 
    if(compressed_)
    {
      // in parallel runs check sequence number of dof manager 
      if( this->grid_.comm().size() > 1 )
      {
        if( sequence_ == this->dofManager_.sequence() ) return false;
      }
      else
      {
        // for serial runs just return 
        return false;
      }
    }

    // mark all leaf elements 
    markAllUsed(); 

    // true if a least one dof must be copied 
    bool haveToCopy = codimLeafSet_[0].compress(); 
    if(higherCodims_)
    {
      for(int i=1; i<ncodim; ++i) 
        haveToCopy = (codimLeafSet_[i].compress()) ? true : haveToCopy;
    }

    // now status is compressed 
    compressed_ = true;
    // update sequence number 
    sequence_ = this->dofManager_.sequence();
    return haveToCopy;
  }

private:
  void assertCodimSetSize ( int codim ) const
  {
#ifndef NDEBUG
    const CodimIndexSetType &codimSet = codimLeafSet_[ codim ];
    if( codimSet.realSize() != hIndexSet_.size( codim ) )
    {
      std::cerr << std::endl;
      std::cerr << "Error in AdaptiveLeafIndexSet: "
                << "Real size of Index set for codim " << codim
                << " is wrong." << std::endl;
      std::cerr << "                               "
                << "Real size of codim index set: " << codimSet.realSize()
                << std::endl;
      std::cerr << "                               "
                << "Size of codim " << codim << " in hierarchic index set: "
                << hIndexSet_.size( codim )
                << std::endl;
      abort();
    }
#endif
  }

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

    assertCodimSetSize( codim );
    const CodimIndexSetType &codimSet = codimLeafSet_[ codim ];
    const int hIdx = hIndexSet_.template index( entity );
    const int idx = codimSet.index( hIdx );
    assert( (idx >= 0) && (idx < codimSet.size()) );
    return idx;
  }

  template< int codim >
  IndexType
  subIndex ( const typename GridType :: template Codim< 0 > :: Entity &entity,
             int subNumber ) const
  {
    if( (codim != 0) && !codimUsed_[ codim ] )
      setUpCodimSet< codim >();

    assertCodimSetSize( codim );
    const CodimIndexSetType &codimSet = codimLeafSet_[ codim ];
    const int hIdx = hIndexSet_.template subIndex< codim >( entity, subNumber );
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
  //! memorise index 
  // --insertIndex
  void insertIndex ( const ElementType &entity )
  {
    const int index = hIndexSet_.index( entity ) ;
    {
#if HAVE_MPI 
      // we need special treatment for ghosts 
      // ghosts should not be inlcuded in holes list 
      if(entity.partitionType() == GhostEntity)
      {
        codimLeafSet_[ 0 ].insertGhost( index );
        if( higherCodims_ )
        {
          const bool skipGhosts = (pitype != All_Partition);
          ForLoop< InsertGhostSubEntities, 1, dimension > :: 
            apply( *this, entity , skipGhosts );
        }
      }
      else 
#endif
      {
        codimLeafSet_[ 0 ].insert( index );
        if( higherCodims_ )
          ForLoop< InsertSubEntities, 1, dimension > :: apply( *this, entity );
      }

      // now compression is not guaranteed 
      compressed_ = false;
    }
    assert( codimLeafSet_[0].exists( index ) );
  }

  //! set indices to unsed so that they are cleaned on compress  
  // --removeIndex
  void removeIndex( const ElementType &entity )
  {
    // remove entities (only mark that they ar not used anymore)
    codimLeafSet_[0].remove ( hIndexSet_.index( entity ) );

    // skip higher codims (will be removed only compression if not
    // existent anymore) 

    // now compressed state is not guaranteed anymore 
    compressed_ = false;
  }

  // insert index if entities lies below used entity, return 
  // false if not , otherwise return true
  bool insertNewIndex ( const ElementType &entity, bool isLeaf, bool canInsert )
  {
    // if entity isLeaf then we insert index 
    if( isLeaf )
    {
      insertIndex( entity );
      return true;
    }
    
    // canInsert is true if we reached an entity 
    // that already has an index 
    if( canInsert ) 
    {
      // we insert to get an index 
      insertIndex( entity );
      // we remove to unmark, because this is not a leaf entity 
      removeIndex( entity );
      return true;
    }
    else 
    {
      // if we have an valid index, then all children may  also appear in the set 
      // from now on, indices can be inserted 
      return codimLeafSet_[0].validIndex( hIndexSet_.index ( entity ) );
    }
  }

  //! mark indices that are still used and give new indices to 
  //! elements that need one 
  void markAllUsed () 
  {
    // make correct size of vectors 
    resizeVectors();

    // unset all indices 
    for(int i=0; i<ncodim; ++i) 
      if(codimUsed_[i]) codimLeafSet_[i].set2Unused(); 
    
    const PartitionIteratorType pt = All_Partition;
    //const PartitionIteratorType pt = pitype;

    typedef typename GridType:: template Codim<0> :: 
      template Partition<pt> :: LeafIterator LeafIteratorType; 

    // walk over leaf level on locate all needed entities  
    LeafIteratorType endit  = this->grid_.template leafend<0,pt>   ();
    for(LeafIteratorType it = this->grid_.template leafbegin<0,pt> (); 
        it != endit ; ++it )
    {
      this->insertIndex( *it );
    }
  }
  
  //! mark indices that are still used and give new indices to 
  //! elements that need one 
  template <int codim>
  void setUpCodimSet () const
  {
    //std::cout << "Setting up codim " << codim << "\n";
    // resize if necessary 
    codimLeafSet_[codim].resize( hIndexSet_.size(codim) );
    
    typedef typename GridType:: template Codim<codim> :: 
      template Partition<pitype> :: LeafIterator LeafIteratorType; 
    // walk over leaf level on locate all needed entities  
    LeafIteratorType endit  = this->grid_.template leafend<codim,pitype>  ();
    for(LeafIteratorType it = 
        this->grid_.template leafbegin<codim,pitype>(); 
        it != endit ; ++it )
    {
      codimLeafSet_[codim].insert( hIndexSet_.index ( *it ) );
    }

    codimUsed_[codim] = true;
    higherCodims_ = true;
  }

  //! give all entities that lie below the old entities new numbers 
  //! here we need the hierarchic iterator because for example for some
  //! grid more the one level of new elements can be created during adaption 
  //! there for we start to give new number for all elements below the old
  //! element 
  void markAllBelowOld () 
  {
    // mark all existing entries as unused 
    for(int i=0; i<ncodim; ++i) 
    {
      if(codimUsed_[i])
      {
        codimLeafSet_[i].set2Unused(); 
      }
    }
    
    const PartitionIteratorType pt = All_Partition;
    //const PartitionIteratorType pt = pitype;

    // get macro iterator 
    typedef typename GridType::
      template Codim<0>::template Partition<pt>:: LevelIterator LevelIteratorType;

    // iterate over macro level and check all entities hierachically 
    const LevelIteratorType macroend = this->grid_.template lend  <0,pt> (0);
    for(LevelIteratorType macroit = this->grid_.template lbegin<0,pt> (0);
        macroit != macroend; ++macroit )
    {
      checkEntity( *macroit , false );
    } // end grid walk trough
  }


  //! check whether entity can be inserted or not 
  void checkEntity(const ElementType& entity, const bool wasNew )
  {
    typedef typename ElementType :: HierarchicIterator HierarchicIteratorType;

    // check whether we can insert or not 
    const bool isNew = insertNewIndex( entity , entity.isLeaf() , wasNew );

    // if entity is not leaf go deeper 
    if( ! entity.isLeaf() )
    {
      const int level = entity.level() + 1;

      // got to next level 
      const HierarchicIteratorType endit  = entity.hend   ( level );
      for(HierarchicIteratorType it = entity.hbegin ( level ); it != endit ; ++it )
      {
        checkEntity( *it , isNew );
      }
    }
  }

  //! count elements by iterating over grid and compare 
  //! entities of given codim with given type 
  template< int codim >
  inline int countElements ( GeometryType type ) const;
  
  // print interal data, for debugging only 
  // print if only done, if DEBUG_LEAFINDEXSET is defined 
  void print (const char * msg, bool oldtoo = false ) const;

public:

  // write indexset to xdr file 
  bool write_xdr(const std::string filename, int timestep) 
  {
    const char *path = "";
    std::string fnstr = genFilename(path,filename, timestep);

    // create write stream 
    XDRWriteStream xdr(fnstr);
    
    // write new verion tag 
    int newVerion = myVersionTag;
    xdr.inout( newVerion );
    
    int type = myType;
    xdr.inout( type );

    bool success = true;

    // write whether codim is used 
    for(int i=0; i<ncodim; ++i) 
    {
      success |= xdr.inout( codimUsed_[i] );
    }

    // write all sets 
    for(int i=0; i<ncodim; ++i) 
      success |= codimLeafSet_[i].processXdr(xdr);
    
    return success;
  }

  //! read index set from given xdr file 
  bool read_xdr(const std::string filename , int timestep)
  {
    const char *path = "";
    std::string fnstr = genFilename(path,filename, timestep);
    
    // create read stream 
    XDRReadStream xdr( fnstr );

    // check type 
    int newVersionTag = myVersionTag;
    xdr.inout( newVersionTag );

    // true if given file has new version 
    const bool newVersion = (newVersionTag == myVersionTag);

    bool success = true;

    // if newVersionTag == myType then older version and 
    // scipt reading of codimUsed 
    int type = (newVersion) ? myType : newVersionTag;

    // if new version the read type, 
    // otherwise newVersionTag is the type info
    if( newVersion )
    {
      // read type 
      xdr.inout( type );
    }

    // index set type check  
    if( (type != 2) && (type != myType) )
    {
      DUNE_THROW(InvalidStateException,"AdaptiveLeafIndexSet::read_xdr: wrong type " << type << " given! Expected " << myType);
    }

    // read codim used 
    if( newVersion )
    {
      // read codim is usage 
      for(int i=0; i<ncodim; ++i) 
      {
        success |= xdr.inout( codimUsed_[i] );
      }
    }
    else 
    {
      // set to used by default 
      for(int i=0; i<ncodim; ++i) 
      {
        codimUsed_[i] = true; 
      }
    }
    
    if(type == 2) 
      success |= codimLeafSet_[0].processXdr(xdr);
    else 
    {
      for(int i=0; i<ncodim; ++i) 
        success |= codimLeafSet_[i].processXdr(xdr);
    }

    // in parallel runs we have to compress here
    if(this->grid_.comm().size() > 1)
    {
      // mark for compress 
      compressed_ = false;
    }
    
    return success;
  }

}; // end of class AdaptiveLeafIndexSet 



template< class GridType, PartitionIteratorType pitype >
template< int codim >
inline int AdaptiveLeafIndexSet< GridType , pitype >
  :: countElements ( GeometryType type ) const
{
  typedef typename GridType :: template Codim< codim >
    :: template Partition< pitype > :: LeafIterator
    IteratorType;

  const GridType &grid = this->grid_;

  int count = 0;
  const IteratorType begin = grid.template leafbegin< codim, pitype >();
  const IteratorType end = grid.template leafend< codim, pitype >();
  for( IteratorType it = begin; it != end; ++it )
  {
    if( it->type() == type )
      ++count;
  }
  return count; 
}

} // end namespace Dune 

#endif
