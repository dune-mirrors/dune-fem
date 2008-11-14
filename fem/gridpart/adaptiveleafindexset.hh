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

  enum INDEXSTATE { NEW, USED, UNUSED };

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
      const HIndexSetType &hIndexSet = indexSet.hIndexSet_;
      if( !indexSet.codimUsed_[ codim ] )
        return;
      
      CodimIndexSetType &codimSet = indexSet.codimLeafSet_[ codim ];
      for( int i = 0; i < entity.template count< codim >(); ++i )
        codimSet.insert( hIndexSet.template subIndex< codim >( entity, i ) );
    }
  };

  template< int codim >
  struct RemoveSubEntities
  {
    static void apply ( ThisType &indexSet, const ElementType &entity )
    {
      const HIndexSetType &hIndexSet = indexSet.hIndexSet_;
      if( !indexSet.codimUsed_[ codim ] )
        return;
      
      CodimIndexSetType &codimSet = indexSet.codimLeafSet_[ codim ];
      for( int i = 0; i < entity.template count< codim >(); ++i )
        codimSet.remove( hIndexSet.template subIndex< codim >( entity, i ) );
    }
  };

  // my type, to be revised 
  enum { myType = 6 };
  enum { myVersionTag = -665 };

  const HIndexSetType &hIndexSet_;
  mutable CodimIndexSetType codimLeafSet_[ ncodim ];
  // flag for codim is in use or not 
  mutable bool codimUsed_ [ncodim];
  
  // true if all entities that we use are marked as USED 
  bool marked_;

  // true if the used entities were marked by grid walkthrough 
  bool markAllU_;

  // true if any of the higher codims is used 
  mutable bool higherCodims_;

  //! flag is tru if set is in compressed status
  mutable bool compressed_;

public:
  //! type traits of this class (see defaultindexsets.hh)
  typedef DefaultLeafIteratorTypes<GridType> Traits; 

  //! Constructor
  AdaptiveLeafIndexSet (const GridType & grid) 
    : BaseType(grid) 
    , hIndexSet_( SelectorType::hierarchicIndexSet(grid) ) 
    , marked_ (false) , markAllU_ (false) 
    , higherCodims_ (false) // higherCodims are not used by default 
    , compressed_(true) // at start the set is compressed 
  {
    // codim 0 is used by default
    codimUsed_[0] = true;

    // all higher codims are not used by default
    for(int i=1; i<ncodim; ++i) codimUsed_[i] = false;
    
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
    marked_ = true;
  }

  //! Unregister entity from set 
  void removeEntity(const typename GridType::template Codim<0>::Entity & en )
  {
    this->removeIndex( en ); 
    marked_ = true;
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
  bool compress ()
  {
    if(compressed_) return false;

    //std::cout << marked_ << " m|u " << markAllU_ << "\n";
    // if not marked, mark which indices are still used 
    //if( (!marked_ && markAllU_) || higherCodims_ ) 
    { 
      //std::cout << "Marking the low level for " << this << "\n";
      markAllUsed(); 
    }

    // true if a least one dof must be copied 
    bool haveToCopy = codimLeafSet_[0].compress(); 
    if(higherCodims_)
    {
      for(int i=1; i<ncodim; ++i) 
        haveToCopy = (codimLeafSet_[i].compress()) ? true : haveToCopy;
    }

    // next turn mark again 
    marked_   = false;
    markAllU_ = false;
    
    compressed_ = true;
    return haveToCopy;
  }

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
    const int hIdx = hIndexSet_.template index( entity );
    const int idx = codimLeafSet_[ codim ].index( hIdx );
    assert( idx >= 0 );
    return idx;
  }

  template< int codim >
  IndexType
  subIndex ( const typename GridType :: template Codim< 0 > :: Entity &entity,
             int subNumber ) const
  {
    if( (codim != 0) && !codimUsed_[ codim ] )
      setUpCodimSet< codim >();
    const int hIdx = hIndexSet_.template subIndex< codim >( entity, subNumber );
    const int idx = codimLeafSet_[ codim ].index( hIdx );
    assert( idx >= 0 );
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
    const int index = hIndexSet_.index( entity );
    if( !codimLeafSet_[ 0 ].exists( index ) )
    {
      codimLeafSet_[ 0 ].insert( index );
      if( higherCodims_ )
        ForLoop< InsertSubEntities, 1, dimension > :: apply( *this, entity );
    }
    compressed_ = false;
  }

  //! set indices to unsed so that they are cleaned on compress  
  // --removeIndex
  void removeIndex( const ElementType &entity )
  {
    // if state is NEW or USED the index of all entities is removed 
    const int index = hIndexSet_.index( entity );
    if( codimLeafSet_[ 0 ].exists( index ) )
    {
      codimLeafSet_[0].remove( index );
      if( higherCodims_ )
        ForLoop< RemoveSubEntities, 1, dimension > :: apply( *this, entity );
    }
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
      // this is the case if we haven't reached an entity which already has a number
      // if index >= 0, then all children may also appear in the set 
      // from now on, indices can be inserted 
      return (codimLeafSet_[ 0 ].index( hIndexSet_.index( entity ) ) >= 0);
    }
  }

  //! mark indices that are still used and give new indices to 
  //! elements that need one 
  void markAllUsed () 
  {
    // unset all indices 
    for(int i=0; i<ncodim; i++) 
      if(codimUsed_[i]) codimLeafSet_[i].set2Unused(); 
    
    typedef typename GridType:: template Codim<0> :: 
      template Partition<pitype> :: LeafIterator LeafIteratorType; 
    // walk over leaf level on locate all needed entities  
    LeafIteratorType endit  = this->grid_.template leafend<0,pitype>   ();
    for(LeafIteratorType it = this->grid_.template leafbegin<0,pitype> (); 
        it != endit ; ++it )
    {
      this->insertIndex( *it );
    }
    marked_ = true;
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

    //codimLeafSet_[codim].print("setup codim");
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
    typedef typename GridType:: template Codim<0> :: 
      template Partition<pitype> :: LevelIterator LevelIteratorType; 

    int maxlevel = this->grid_.maxLevel();
   
    for(int i=0; i<ncodim; i++) 
    {
      if(codimUsed_[i])
      {
        codimLeafSet_[i].set2Unused(); 
      }
    }
    
    for(int level = 0; level<=maxlevel; level++)
    {
      LevelIteratorType levelend    = 
        this->grid_.template lend  <0,pitype> (level);
      for(LevelIteratorType levelit = 
          this->grid_.template lbegin<0,pitype> (level);
          levelit != levelend; ++levelit )
      {
        typedef typename GridType::template Codim<0>::
              Entity::HierarchicIterator HierarchicIteratorType; 
       
        // if we have index all entities below need new numbers 
        bool areNew = false; 

        // check whether we can insert or not 
        areNew = insertNewIndex ( *levelit , levelit->isLeaf() , areNew ); 
        
        HierarchicIteratorType endit  = levelit->hend   ( level + 1 );
        for(HierarchicIteratorType it = levelit->hbegin ( level + 1 ); it != endit ; ++it )
        {
          // areNew == true, then index is inserted 
          areNew = insertNewIndex  ( *it , it->isLeaf() , areNew ); 
        }

      } // end grid walk trough
    } // end for all levels 

    // means on compress we have to mark the leaf level 
    marked_ = false;
    markAllU_ = true;
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
