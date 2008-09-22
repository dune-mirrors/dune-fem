#ifndef DUNE_IDBASEDLEAFINDEXSET_HH
#define DUNE_IDBASEDLEAFINDEXSET_HH

//- system includes 
#include <map>
#include <stack>

//- local includes 
#include "adaptiveleafindexset.hh"

namespace Dune { 

//! Send vector to output stream
template<class MapType>
void printMap (std::ostream& s, const MapType& m)
{
  typedef typename MapType :: const_iterator iterator;
  iterator end = m.end();
  for(iterator it=m.begin(); it != end; ++it) 
  {
    int idx = (*it).second;
    s << idx << " ";  
  }
  s << std::endl;
}
  
//! \brief LeafIndexSet based on local ids using std::map to build mapping
//!from id to index 
template <class IdSetImp, class IndexSetImp, int codim>  
class IdBasedCodimIndexSet
{
  typedef IdSetImp IdSetType;
  typedef IndexSetImp IndexSetType;
  typedef typename IdSetType :: IdType IdType;
  
private:
  typedef std::map< IdType, int  > IndexStorageType;
  typedef std::map< int , IdType > IdStorageType;
  typedef std::vector < int > HolesIndicesType;

  typedef MutableArray<int> IndexVectorType;
  
  // the mapping of the global to leaf index 
  mutable IndexStorageType leafIndex_;
  mutable IdStorageType leafId_;

  // old indices 
  mutable IndexStorageType oldLeafIndex_;

  // list of old and new index for all holes 
  mutable IndexVectorType oldIndexVec_;
  mutable IndexVectorType newIndexVec_;

  const IdSetType    & idSet_;
  const IndexSetType & indexSet_;
 
  // next index to give away 
  int nextFreeIndex_;

  int oldIds_;

  // no copying 
  IdBasedCodimIndexSet (const IdBasedCodimIndexSet& );

public:
  //! Constructor
  template< class Iterator >
  IdBasedCodimIndexSet ( const IdSetType &idSet,
                         const IndexSetType &indexSet,
                         const Iterator &begin,
                         const Iterator &end )
  : idSet_( idSet ),
    indexSet_( indexSet ),
    nextFreeIndex_( indexSet_.size( codim ) )
  {
    createMaps( begin, end );

    // set memory over estimation 
    oldIndexVec_.setMemoryFactor( 1.1 );
    newIndexVec_.setMemoryFactor( 1.1 );
  }

  template< class Iterator >
  void createMaps ( const Iterator &begin, const Iterator &end );

  //! make to index numbers consecutive 
  //! return true, if at least one hole was closed 
  void resize ()
  {
    nextFreeIndex_ = indexSet_.size(codim);
    oldLeafIndex_.clear();
    oldLeafIndex_ = leafIndex_;

    leafId_.clear();
    leafIndex_.clear();
    oldIds_ = 0;
  }

  //! make to index numbers consecutive 
  //! return true, if at least one hole was closed 
  bool compress ()
  {
    int actSize = leafIndex_.size ();
    nextFreeIndex_ = actSize;

    // create holes index vector 
    HolesIndicesType holesIdx(oldLeafIndex_.size(), -1);   

    int noHole = 0;
    {
      // remove all indices that are bigger than nextFreeIndex 
      typedef typename IndexStorageType :: const_iterator iterator;
      iterator end = oldLeafIndex_.end();
      for(iterator it = oldLeafIndex_.begin(); it != end; ++it) 
      {
        int idx = (*it).second;
        assert( idx >= 0 );
        
        // check if index is valid and is a hole 
        if( (idx < actSize) && (idx >= 0) )
        {
          // if id is not in the new set, then its a hole 
          if( leafId_.find(idx) == leafId_.end() )
          {
            holesIdx[noHole] = idx;
            ++noHole;
          }
        }
      }
    }

    assert( (leafIndex_.size() != (size_t) actSize) ? 
        (std::cerr << actSize << " s|ls " << leafIndex_.size() << "\n",0) : 1);

    {
      // resize hole lists 
      oldIndexVec_.resize( actSize );
      newIndexVec_.resize( actSize );
      
      int hole = 0;
      typedef typename IndexStorageType :: iterator iterator;
      const iterator end = leafIndex_.end();
      for(iterator it = leafIndex_.begin(); it != end; ++it) 
      {
        const int idx = (*it).second;
        if( idx >= nextFreeIndex_ ) 
        {
          assert( hole < oldIndexVec_.size() );

          // put to old idx list 
          oldIndexVec_[hole] = idx;
          
          int newIdx = -1; 
          if( hole >= noHole )
          {
            newIdx = nextFreeIndex_; 
            ++nextFreeIndex_; 
          }
          else 
          {
            newIdx = holesIdx[hole];
          }
          
          assert( hole < newIndexVec_.size() );
          
          // store new index 
          newIndexVec_[hole] = newIdx;
          
          // set new index to index map 
          (*it).second = newIdx;

          leafId_.erase(idx); 
          leafId_[newIdx] = (*it).first;

          // next hole 
          ++hole;
        }

        assert( (*it).second < size() ); 
      }

      // set current number of holes
      oldIndexVec_.resize( hole );
      newIndexVec_.resize( hole );
    }

    // clear old values 
    oldLeafIndex_.clear();
    
    // check that index set is consecutive 
    // only done in debug mode 
    checkConsecutive();
    return true; 
  }

  void checkConsecutive()
  {
#ifndef NDEBUG
    typedef typename IndexStorageType :: const_iterator iterator;
    iterator end    = leafIndex_.end();
    for(iterator it = leafIndex_.begin(); it != end; ++it) 
    {
      int idx = (*it).second;
      assert( idx < nextFreeIndex_ );
    }
#endif
  }

  //! return how much extra memory is needed for restriction 
  int additionalSizeEstimate () const { return size(); }

  //! return size of grid entities per level and codim 
  int size () const
  {
    return nextFreeIndex_;
  }
  
  //! return size of grid entities per level and codim 
  int realSize () const
  {
    return size();
  }

  //! return leaf index for given hierarchic number  
  template <class EntityType> 
  int index ( const EntityType & en, int ) const
  {
    // assert if index was not set yet 
    assert( exists( en ) );
    assert( leafIndex_ [ idSet_.id(en) ] < size() );
    return leafIndex_ [ idSet_.id(en) ];
  }
 
  //! return state of index for given hierarchic number  
  template <class EntityType> 
  bool exists ( const EntityType & en) const
  {
    return (leafIndex_.find( idSet_.id(en) ) != leafIndex_.end()); 
  }
 
  //! return number of existing holes 
  int numberOfHoles () const
  {
    assert( oldIndexVec_.size() == newIndexVec_.size() );
    return oldIndexVec_.size();
  }

  //! return old index, for dof manager only 
  int oldIndex ( const int idx ) const
  {
    return oldIndexVec_[ idx ];
  }

  //! return new index, for dof manager only returns index 
  int newIndex ( const int idx ) const
  {
    return newIndexVec_[ idx ];
  }
  
  // check index   
  template <class EntityType> 
  void checkIndex (const EntityType & en) 
  {
    IdType id = idSet_.id(en);
    if(leafIndex_.find(id) == leafIndex_.end())
    {
      if(oldLeafIndex_.find(id) == oldLeafIndex_.end())
      {
        int idx = nextFreeIndex_;
        ++nextFreeIndex_;
        leafIndex_[id] = idx; 
        leafId_[idx] = id;
      }
      else 
      {
        ++oldIds_;
        int idx = oldLeafIndex_[id];
        leafIndex_[id]  = idx; 
        leafId_   [idx] = id;
      }
    }
  }
  
  // insert element and create index for element number  
  template <class EntityType> 
  void insert (const EntityType & en ) 
  {
    IdType id = idSet_.id(en);
    if(leafIndex_.find(id) == leafIndex_.end())
    {
      leafIndex_[id] = nextFreeIndex_;
      leafId_[nextFreeIndex_] = id;
      ++nextFreeIndex_;
    }
  }
  
  // remove index actually is done in compression 
  template <class EntityType> 
  void remove ( const EntityType & en ) 
  {
  }
  
  void print (const char * msg, bool oldtoo = false ) const {};
  
  // read/write from/to xdr stream 
  bool processXdr(XDR *xdrs)
  {
    /*
    xdr_int ( xdrs, &nextFreeIndex_ );
    xdr_int ( xdrs, &actSize_ );
    leafIndex_.processXdr(xdrs);
    state_.processXdr(xdrs);
    return true;
    */
    return false;
  }
}; // end of class 


  template< class IdSetImp, class IndexSetImp, int codim >
  template< class Iterator >
  void IdBasedCodimIndexSet< IdSetImp, IndexSetImp, codim >
    :: createMaps ( const Iterator &begin, const Iterator &end )
  {
    for( Iterator it = begin; it != end; ++it )
    {
      IdType id = idSet_.id(*it);
      if(leafIndex_.find(id) == leafIndex_.end())
      {
        int idx = indexSet_.index(*it);
        leafIndex_[id] = idx;    
        leafId_  [idx] = id;
      }
    }

    checkConsecutive();
  }



/** \class IdBasedLeafIndexSet
 *  \brief consecutive leaf index set based on the grid's local id set
 *
 *  This index set maps an entity to a consecutive leaf index. Internally,
 *  the entity's local id is mapped to the leaf index using an STL map.
 *
 *  \note Only codim 0 is working at the moment
 *
 *  \note Local id's are only unique within the InteriorBorder_Partition. Dune
 *        requires the ids of periodic ghost and overlap entities to equal the
 *        ids of the corresponding interior element (the index has to be unique,
 *        though). So be careful when using this index set with a periodic grid.
 *
 *  \todo Make the index set work for all codimensions
 */
template< class Grid >
class IdBasedLeafIndexSet
: public ConsecutivePersistentIndexSet
  < Grid, IdBasedLeafIndexSet< Grid >, DefaultLeafIteratorTypes< Grid > >
{
  typedef IdBasedLeafIndexSet< Grid > ThisType;
  typedef ConsecutivePersistentIndexSet
    < Grid, ThisType, DefaultLeafIteratorTypes< Grid > >
    BaseType;

public:
  typedef Grid GridType;

  static const int dim = GridType :: dimension;

  //! type of index 
  typedef typename BaseType :: IndexType IndexType;

private:
  friend class Conversion< ThisType, EmptyIndexSet >;
  
  // my type, to be revised 
  enum { myType = 8 };

  typedef typename GridType :: Traits :: LeafIndexSet LeafIndexSetType;
  typedef typename GridType :: Traits :: LocalIdSet LocalIdSetType;
  typedef typename GridType :: template Codim<0> :: Entity EntityCodim0Type;

  typedef IdBasedCodimIndexSet< LocalIdSetType, LeafIndexSetType, 0 >
    IdBasedCodimIndexSetType;

  const LeafIndexSetType &leafSet_; 
  const LocalIdSetType &idSet_; 


  mutable IdBasedCodimIndexSetType pLeafSet_;
  
  // true if all entities that we use are marked as USED 
  bool compressed_;

  // no copying, no assignment
  IdBasedLeafIndexSet ( const ThisType & );
  ThisType &operator= ( const ThisType & );

public:
  //! type traits of this class
  typedef DefaultLeafIteratorTypes<GridType> Traits; 

  //! Constructor
  inline explicit IdBasedLeafIndexSet ( const GridType &grid )
  : BaseType( grid ),
    leafSet_( grid.leafIndexSet() ),
    idSet_( grid.localIdSet() ),
    pLeafSet_( idSet_, leafSet_,
               grid.template leafbegin< 0 >(), grid.template leafend< 0 >() ),
    compressed_( true )
  {}

  inline explicit IdBasedLeafIndexSet ( const GridType *const grid )
  : BaseType( *grid ),
    leafSet_( grid->leafIndexSet() ),
    idSet_( grid->localIdSet() ),
    pLeafSet_( idSet_, leafSet_,
               grid->template leafbegin< 0 >(), grid->template leafend< 0 >() ),
    compressed_( true )
  {}

  //! Destructor
  virtual ~IdBasedLeafIndexSet () {};

  //! return type of index set, for GrapeDataIO
  int type () const { return myType; }

  //! return name of index set, for GrapeDataIO
  std::string name () const { return "IdBasedLeafIndexSet"; }

  //****************************************************************
  //
  //  INTERFACE METHODS for DUNE INDEX SETS 
  //
  //****************************************************************
  //! return size of grid entities per level and codim 
  int size ( int codim , GeometryType type ) const
  {
    return persistentLeafSet().size();
  }

  //! \brief return size of grid entities per geometry type
  int size (GeometryType type) const
  {
    // return size using the other method
    return size(dim - type.dim(), type);
  }

  //! return size of grid entities per level and codim 
  int size ( int codim ) const
  {
    size_t s = 0;
    const std::vector< GeometryType > & geomTypes = leafSet_.geomTypes(codim);
    for(size_t i=0; i<geomTypes.size(); ++i)
      s += size(codim,geomTypes[i]);
    return s;
  }
  
  //! returns vector with geometry tpyes this index set has indices for
  const std::vector <GeometryType> & geomTypes (int codim) const 
  {
    return leafSet_.geomTypes(codim);
  }
 
#ifdef INDEXSET_HAS_ITERATORS
  /** @brief Iterator to one past the last entity of given codim for partition type
   *  Here the grids leaf iterator is used 
   */
  template<int cd, PartitionIteratorType pitype>
  typename DefaultLeafIteratorTypes<GridType>::template Codim<cd>::
    template Partition<pitype>::Iterator end () const
  {
    return leafSet_.template end<cd,pitype> ();
  }

  /** @brief Iterator to first entity of given codimension and partition type.
   *  Here the grids leaf iterator is used 
   */
  template<int cd, PartitionIteratorType pitype>
  typename DefaultLeafIteratorTypes<GridType>::template Codim<cd>::
    template Partition<pitype>::Iterator begin () const
  {
    return leafSet_.template begin<cd,pitype> ();
  }
#endif

  //! \brief returns true if entity is contained in index set 
  template <class EntityType>
  bool contains (const EntityType & en) const
  {
    return persistentLeafSet().exists(en); 
  }

  //****************************************************************
  //
  //  METHODS for Adaptation with DofManger 
  //
  //****************************************************************
  //! insert entity into set 
  void insertEntity(const typename GridType::template Codim<0>::Entity & en )  
  {
    this->insertIndex( en );
  }

  //! Unregister entity which will be removed from the grid
  void removeEntity(const typename GridType::template Codim<0>::Entity & en )
  {
    this->removeIndex( en ); 
  }

  void resize () 
  {
    // give all entities that lie below 
    // the old entities new numbers 
    markAllBelowOld ();
  }

  //! make to index numbers consecutive 
  //! return true, if at least one hole was closed 
  bool compress ()
  {
    // if already compressed, do nothing 
    if ( compressed_ ) return false;

    // prepare for resize 
    persistentLeafSet().resize();

    // if not marked, mark which indices are still used 
    markAllUsed(); 

    // true if a least one dof must be copied 
    bool haveToCopy = persistentLeafSet().compress(); 

    // now index set is in compressed state  
    compressed_ = true;
    return haveToCopy;
  }

  //! return approximate size that is used during restriction 
  int additionalSizeEstimate () const 
  { 
    return persistentLeafSet().additionalSizeEstimate();
  }

  //! return global index 
  //! for dof mapper 
  // --index 
  template <int codim, class EntityType>
  IndexType indexImp (const EntityType & en, int num) const
  {
    return persistentLeafSet().index(en,num);
  }
  
  //! return number of exisiting holes 
  int numberOfHoles(const int codim)  const
  {
    return persistentLeafSet().numberOfHoles();
  }

  //! return old index, for dof manager only 
  int oldIndex (const int hole, const int codim ) const
  {
    return persistentLeafSet().oldIndex( hole );
  }

  //! return new index, for dof manager only returns index 
  int newIndex (const int hole , const int codim ) const
  {
    return persistentLeafSet().newIndex( hole );
  }

protected:
  //! memorise index 
  // --insertIndex
  void insertIndex(const EntityCodim0Type & en)
  {
    persistentLeafSet().insert(en);
    compressed_ = false;
  }

  //! set indices to unsed so that they are cleaned on compress  
  // --removeIndex
  void removeIndex(const EntityCodim0Type & en)
  {
    persistentLeafSet().remove(en);
    compressed_ = false;
  }

  // insert index if entities lies below used entity, return 
  // false if not , otherwise return true
  bool insertNewIndex (const EntityCodim0Type & en, bool isLeaf , bool canInsert )
  {
    // if entity isLeaf then we insert index 
    if(isLeaf)
    {
      this->insertIndex(en );
      return true;
    }
    
    // which is the case if we havent reached a entity which has 
    // already a number 
    if(!canInsert) 
    {
      // from now on, indices can be inserted 
      if( persistentLeafSet().exists(en) )
      {
        return true; 
      }

      // we have to go deeper 
      return false;
    }
    else 
    {
      this->removeIndex( en );
      // set unused here, because index is only needed for prolongation 
    }
    return true;
  }

  IdBasedCodimIndexSetType & persistentLeafSet() 
  { 
    return pLeafSet_; 
  }
  const IdBasedCodimIndexSetType & persistentLeafSet() const 
  { 
    return pLeafSet_; 
  }

  //! mark indices that are still used and give new indices to 
  //! elements that need one 
  void markAllUsed () 
  {
    // unset all indices 
    typedef typename GridType:: template Codim<0> :: LeafIterator LeafIteratorType; 
    // walk over leaf level on locate all needed entities  
    LeafIteratorType endit  = this->grid_.template leafend<0>   ();
    for(LeafIteratorType it = this->grid_.template leafbegin<0> (); 
        it != endit ; ++it )
    {
      persistentLeafSet().checkIndex( *it );
    }
  }
  
  //! give all entities that lie below the old entities new numbers 
  //! here we need the hierarchic iterator because for example for some
  //! grid more the one level of new elements can be created during adaption 
  //! there for we start to give new number for all elements below the old
  //! element 
  void markAllBelowOld () 
  {
    typedef typename GridType::
      template Codim<0>::LevelIterator LevelIteratorType; 

    // iterate over macro level and check all entities hierachically 
    const LevelIteratorType macroend = this->grid_.template lend  <0> (0);
    for(LevelIteratorType macroit = this->grid_.template lbegin<0> (0);
        macroit != macroend; ++macroit )
    {
      checkEntity( *macroit , false );
    } // end grid walk trough
  }

  //! check whether entity can be inserted or not 
  void checkEntity(const EntityCodim0Type& en, const bool wasNew )
  {
    typedef typename EntityCodim0Type :: HierarchicIterator HierarchicIteratorType;

    // check whether we can insert or not 
    const bool isNew = insertNewIndex ( en , en.isLeaf() , wasNew );

    // if entity is not leaf go deeper 
    if( ! en.isLeaf() )
    {
      const int level = en.level() + 1;

      // got to next level 
      const HierarchicIteratorType endit  = en.hend   ( level );
      for(HierarchicIteratorType it = en.hbegin ( level ); it != endit ; ++it )
      {
        checkEntity( *it , isNew );
      }
    }
  }


  //! count elements by iterating over grid and compare 
  //! entities of given codim with given type 
  template< int codim >
  int countElements ( GeometryType type ) const;
  
  // print interal data, for debugging only 
  // print if only done, if DEBUG_LEAFINDEXSET is defined 
  void print (const char * msg, bool oldtoo = false ) const;

public:

  // write indexset to xdr file 
  bool write_xdr(const std::basic_string<char> filename, int timestep) 
  {
    /*
    assert( false );
    abort();
    FILE  *file;
    XDR   xdrs;
    const char *path = "";

    std::basic_string<char> fnstr = genFilename(path,filename, timestep);
    const char * fn = fnstr.c_str();
    file = fopen(fn, "wb");
    if (!file)
    {
      std::cerr << "\aERROR in AdaptiveLeafIndexSet::write_xdr(..): couldnot open " << filename << std::endl;
      std::cerr.flush();
      return false;
    }

    xdrstdio_create(&xdrs, file, XDR_ENCODE);
    int type = myType;
    xdr_int ( &xdrs, &type );
    if(type != myType)
    {
      std::cerr << "\nERROR: AdaptiveLeafIndexSet: wrong type choosen! \n\n";
      assert(type == myType);
    }

    for(int i=0; i<=dim; i++) codimLeafSet_[i].processXdr(&xdrs);

    xdr_destroy(&xdrs);
    fclose(file);
    */
    return true;
  }

  //! read index set from given xdr file 
  bool read_xdr(const std::basic_string<char> filename , int timestep)
  {
    assert( false );
    abort();
    /*
    FILE   *file;
    XDR     xdrs;
    const char *path = "";

    std::basic_string<char> fnstr = genFilename(path,filename, timestep);
    const char * fn = fnstr.c_str();
    std::cout << "Reading <" << fn << "> \n";
    file = fopen(fn, "rb");
    if(!file)
    {
      std::cerr <<"\aERROR in AdaptiveLeafIndexSet::read_xdr(..): couldnot open <%s>!\n" << filename << std::endl;
      std::cerr.flush();
      return(false);
    }

    // read xdr 
    xdrstdio_create(&xdrs, file, XDR_DECODE);
    
    int type = myType;
    xdr_int ( &xdrs, &type );
    if( (type != 2) && (type != myType) )
    {
      std::cerr << "\nERROR: AdaptiveLeafIndexSet: wrong type choosen! \n\n";
      assert(type == myType);
    }

    if(type == 2) 
      codimLeafSet_[0].processXdr(&xdrs);
    else 
    {
      for(int i=0; i<=dim; i++) codimLeafSet_[i].processXdr(&xdrs);
    }

    xdr_destroy(&xdrs);
    fclose(file);
    */

    print("read Index set ");
    return true;
  }

}; // end of class AdaptiveLeafIndexSet 



template< class GridType >
template< int codim >
inline int IdBasedLeafIndexSet< GridType >
  :: countElements ( GeometryType type ) const
{
  typedef typename GridType :: template Codim< codim >
    :: template Partition< All_Partition > :: Iterator
    IteratorType;

  const GridType &grid = this->grid_;

  int count = 0;
  const IteratorType begin = grid.template leafbegin< codim, All_Partition >();
  const IteratorType end = grid.template leafend< codim, All_Partition >();
  for( IteratorType it = begin; it != end; ++it )
  {
    if( it->type() == type )
      ++count;
  }
  return count; 
}



template <class GridType> 
inline void IdBasedLeafIndexSet<GridType>::
print (const char * msg, bool oldtoo ) const 
{
#ifdef DEBUG_LEAFINDEXSET
    const CodimLeafIndexSet & cls = codimLeafSet_[0];
    std::cout << "Size " << cls.size() << "\n";
    std::cout << "i    |   val    | state  \n";
    int actSize =0;
    
    for(int i=0; i< cls.realSize(); i++)
    {
      if( cls.state( i ) != UNUSED ) actSize ++;
      std::cout << i << " | " << cls.index(i) << " | " << cls.state( i );
      std::cout << "\n";
    }
    
    std::cout << "Real Size " << cls.size() << "\n";
    std::cout << "ActSize   " << actSize << "\n";
    std::cout << "Grid global Size " << hIndexSet_.size(0) << "\n";

    std::cout << msg;
#endif
}

} // end namespace Dune 

#endif
