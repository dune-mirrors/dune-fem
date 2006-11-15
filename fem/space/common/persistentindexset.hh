#ifndef DUNE_PERSISTENTLEAFINDEXSET_HH
#define DUNE_PERSISTENTLEAFINDEXSET_HH

//- system includes 
#include <map>
#include <stack>

//- local includes 
#include "adaptiveleafindexset.hh"

namespace Dune { 


//***********************************************************************
//
//  Index Set for one codimension
//  --CodimLeafIndexSet 
//
//***********************************************************************

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
  

template <class IdSetImp, class IndexSetImp, int codim>  
class PersistentLeafIndexSet
{
  typedef IdSetImp IdSetType;
  typedef IndexSetImp IndexSetType;
  typedef typename IdSetType :: IdType IdType;
  
private:
  typedef std::map< IdType, int  > IndexStorageType;
  typedef std::map< int , IdType > IdStorageType;
  typedef std::vector < int > HolesIndicesType;
  
  // the mapping of the global to leaf index 
  mutable IndexStorageType leafIndex_;

  mutable IdStorageType leafId_;
  mutable IdStorageType holesId_;

  // old indices 
  mutable IndexStorageType holes_;
  mutable IndexStorageType oldLeafIndex_;

  const IdSetType    & idSet_;
  const IndexSetType & indexSet_;
 
  // next index to give away 
  int nextFreeIndex_;

  int oldIds_;
public:
  //! Constructor
  PersistentLeafIndexSet (const IdSetType & idSet, const IndexSetType & indexSet) 
    : idSet_(idSet), indexSet_(indexSet) , 
      nextFreeIndex_ (indexSet_.size(codim)) 
  {
    createMaps();
  }

  void createMaps () 
  {
    typedef typename IndexSetType :: template Codim<codim> ::
      template Partition<All_Partition> :: Iterator IteratorType; 

    IteratorType end = indexSet_.template end<codim,All_Partition> (); 
    for(IteratorType it = indexSet_.template begin<codim,All_Partition>(); 
        it != end; ++it )
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
    //std::cout << "Created maps \n";
  }

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

    holes_.clear();
    holesId_.clear();
    HolesIndicesType holesIdx(oldLeafIndex_.size());   
    for(size_t i=0; i<holesIdx.size(); ++i)
      holesIdx[i] = -1;

    //std::cout << "Start compressing, found " << oldIds_ << " elements \n";
    //std::cout << "Size of leafindex set = " << indexSet_.size(0) << "\n";

    //assert( actSize == indexSet_.size(codim) );

    int noHole = 0;
    {
      // remove all indices that are bigger than nextFreeIndex 
      typedef typename IndexStorageType :: const_iterator iterator;
      iterator end = oldLeafIndex_.end();
      for(iterator it = oldLeafIndex_.begin(); it != end; ++it) 
      {
        int idx = (*it).second;
        assert( idx >= 0 );
        // normaly one would check >= 0, but does not work for UGGrid
        if( (idx < actSize) && (idx > 0) )
        {
          // if id is not in the new set, then its a hole 
          if( leafId_.find(idx) == leafId_.end() )
          {
            //std::cout << "Found hole " << idx << " \n";
            holesIdx[noHole] = idx;
            ++noHole;
          }
        }
      }
    }

    //std::cout << "Found " << noHole << " holes \n";

    assert( (leafIndex_.size() != (size_t) actSize) ? 
        (std::cerr << actSize << " s|ls " << leafIndex_.size() << "\n",0) : 1);
    {
      int hole = 0;
      typedef typename IndexStorageType :: iterator iterator;
      iterator end = leafIndex_.end();
      for(iterator it = leafIndex_.begin(); it != end; ++it) 
      {
        int idx = (*it).second;
        //std::cout << "Check leaf index " << idx << " \n";
        if( idx >= nextFreeIndex_ ) 
        {
          IdType id = (*it).first; 
          holesId_[hole] = id;
          holes_[id] = idx; 

          int newIdx = -1; 
          if( hole >= noHole )
          {
            //std::cout << "Run out of holes \n";
            newIdx = nextFreeIndex_; 
            ++nextFreeIndex_; 
          }
          else 
          {
            newIdx = holesIdx[hole];
            //std::cout << "got hole " << newIdx << "\n";
          }
          ++hole;
          
          (*it).second = newIdx;
          //leafIndex_.erase(id);
          //leafIndex_[id] = newIdx;

          leafId_.erase(idx); 
          leafId_[newIdx] = id;
        }
        assert( (*it).second < size() ); 
      }
    }

    //std::cout << "We have " << numberOfHoles() << " holes \n";
    oldLeafIndex_.clear();
 
    /*
    {
      typedef typename IndexStorageType :: const_iterator iterator;
      iterator end = leafIndex_.end();
      for(iterator it = leafIndex_.begin(); it != end; ++it) 
      {
        IdType id = (*it).first;
        int idx = (*it).second;
        iterator end2 = leafIndex_.end();
        for(iterator it2 = leafIndex_.begin(); it2 != end; ++it2) 
        {
          IdType id2 = (*it2).first;
          if( id2 != id ) 
          {
            if(idx == (*it2).second)
            {
              std::cout << id << " id and " << id2 << " have index " << idx << "\n";
              assert( idx != (*it2).second); 
            }
          }
          
        }
      }
      
    }
    */
    
    //checkConsecutive();
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
      //if( idx >= nextFreeIndex_ ) 
      //  std::cout << "index to big " << idx << " " << nextFreeIndex_ << "\n";
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
  int index ( const EntityType & en, int num ) const
  {
    // assert if index was not set yet 
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
    return holes_.size();
  }

  //! return old index, for dof manager only 
  int oldIndex ( int idx ) const
  {
    assert( holesId_.find(idx) != holesId_.end() );
    IdType id = holesId_[idx];
    assert( holes_.find(id) != holes_.end() );
    return holes_[id]; 
  }

  //! return new index, for dof manager only returns index 
  int newIndex ( int idx ) const
  {
    assert( holesId_.find(idx) != holesId_.end() );
    IdType id = holesId_[idx];
    assert( leafIndex_.find(id) != leafIndex_.end() ); 
    assert( leafIndex_[id] < size() );
    return leafIndex_[id]; 
  }

  // check index   
  template <class EntityType> 
  void checkIndex (const EntityType & en) 
  {
    IdType id = idSet_.id(en);
    //std::cout << id << " check id \n";
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
        //std::cout << "Found old index " << idx << "\n";
        leafIndex_[id]  = idx; 
        leafId_   [idx] = id;
      }
    }

    //std::cout << "Index of id " << id << " is " << leafIndex_[id] << "\n";
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
  
  // insert element and create index for element number  
  template <class EntityType> 
  void remove ( const EntityType & en ) 
  {
    /*
    IdType id = idSet_.id(en);
    if(leafIndex_.find(id) != leafIndex_.end())
    {
      int idx = leafIndex_[id];
      leafIndex_.erase(id);
      leafId_.erase(idx);  
    }
    */
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

//******************************************************************
//
// Indexset that provides consecutive indicies for the leaf level
// this index set uses the grid hierarchical index  
//
//******************************************************************
/*! 
  This index set generates a consecutive leaf index out of the unique
  global index of each entity. This index set can be used instead of the
  default grid index sets and can be generated for each grid implementation.

  Note that only codim = 0 is working at the moment. 

  Future work. Index set for each codim.
*/


template <class GridType>
class DefaultAdaptiveLeafIndexSet : 
  public IndexSet<GridType, DefaultAdaptiveLeafIndexSet<GridType>, LeafIteratorTypes<GridType> >,
  public DefaultGridIndexSetBase <GridType>
{
public:
  enum { ncodim = GridType::dimension + 1 };

private:

  // busines as usual 


  // count elements of set by iterating the grid 
  template <class AdLeafSet, int codim >
  struct CountElements
  {
    static inline int count (const AdLeafSet & ls , int cd, GeometryType type )
    {
      if( cd == codim ) 
      {
        return ls.template countElements<codim> (type);
      }
      else 
        return CountElements < AdLeafSet, codim-1> :: count (ls,cd,type);
    }
  };

  // count elements of set by iterating the grid 
  template <class AdLeafSet>
  struct CountElements<AdLeafSet,0>
  {
    static inline int count (const AdLeafSet & ls , int cd, GeometryType type )
    {
      enum { codim = 0 };
      if( cd == codim ) 
      {
        return ls.template countElements<codim> (type);
      }
      else 
        return 0;
    }
  };

  //***************************************************************************//
  
  // direct index return from method index (const EntityType & en )
  template <class AdLeafSet, class HSetImp, class EntityType, int enCodim >
  struct DirectIndexWrapper
  {
    static inline int index (const AdLeafSet & ls , const HSetImp & hset, 
                             const EntityType & en, bool cdUsed )
    {
      return hset.index(en);
    }
  };

  //**************************************************************//

  // index return from method index (const EntityType & en, int num )
  // this puts the index method and the subIndex methods together 
  template <class AdLeafSet, class HSetImp, class EntityType, int enCodim, int codim >
  struct IndexWrapper
  {
    static inline int index (const AdLeafSet & ls , const HSetImp & hset, 
                             const EntityType & en , int num )
    {
      return hset.index(en);
    }
  };

  //! if codim > codim of entity use subIndex 
  template <class AdLeafSet, class HSetImp, class EntityType>
  struct IndexWrapper<AdLeafSet,HSetImp,EntityType,0,1>
  {
    static inline int index (const AdLeafSet & ls , const HSetImp & hset, 
                             const EntityType & en , int num )
    {
      enum { codim = 1 };
      return hset.template subIndex<codim> (en , num ) ;
    }
  };

  //! if codim > codim of entity use subIndex 
  template <class AdLeafSet, class HSetImp, class EntityType>
  struct IndexWrapper<AdLeafSet,HSetImp,EntityType,0,2>
  {
    static inline int index (const AdLeafSet & ls , const HSetImp & hset, 
                             const EntityType & en , int num )
    {
      enum { codim = 2 };
      return hset.template subIndex<codim> (en , num ) ;
    }
  };

  //! if codim > codim of entity use subIndex 
  template <class AdLeafSet, class HSetImp, class EntityType>
  struct IndexWrapper<AdLeafSet,HSetImp,EntityType,0,3>
  {
    static inline int index (const AdLeafSet & ls , const HSetImp & hset, 
                             const EntityType & en , int num )
    {
      enum { codim = 3 };
      return hset.template subIndex<codim> (en , num ) ;
    }
  };

  //******************************************************************
  //  partial specialisation for the insertion of all sub entity indices 
  //******************************************************************
  template <class HSetImp, class CodimLeafSet, class EntityType, int codim> 
  struct PartialSpec 
  {
    // only for higher codims 
    CompileTimeChecker< (codim > 1) ? true : false> check; 
    
    static inline void iterateCodims (const HSetImp & hIndexSet , 
        CodimLeafSet (&cls)[ncodim], const EntityType & en , bool (&cdUsed)[ncodim])
    {
      CodimLeafSet & lset = cls[codim];

      // if codim is used then insert all sub entities of this codim 
      if(cdUsed[codim])
      {
        for(int i=0; i<en.template count<codim> (); i++)
        {
          lset.insert( hIndexSet. template subIndex<codim> (en,i) );   
        }
      }
      
      PartialSpec<HSetImp,CodimLeafSet,EntityType,codim-1> :: 
        iterateCodims (hIndexSet, cls, en , cdUsed );
    }
    
    static inline void removeCodims (const HSetImp & hIndexSet , 
        CodimLeafSet (&cls)[ncodim], const EntityType & en , bool (&cdUsed)[ncodim])
    {
      CodimLeafSet & lset = cls[codim];

      // if codim is already used, then also remove entities of this codim
      if(cdUsed[codim])
      {
        for(int i=0; i<en.template count<codim> (); i++)
        {
          lset.remove( hIndexSet. template subIndex<codim> (en,i) );   
        }
      }
      
      PartialSpec<HSetImp,CodimLeafSet,EntityType,codim-1> :: 
        removeCodims (hIndexSet, cls, en , cdUsed );
    }
  };
 
  // specialisation for codim 1 is then end of the loop
  template <class HSetImp, class CodimLeafSet, class EntityType> 
  struct PartialSpec<HSetImp,CodimLeafSet,EntityType,1> 
  {
    static inline void iterateCodims (const HSetImp & hIndexSet , 
        CodimLeafSet (&cls)[ncodim], const EntityType & en , bool (&cdUsed)[ncodim])
    {
      enum { codim = 1 };
      CodimLeafSet & lset = cls[codim];
      // if codim is already used, then also insert entities of this codim
      if(cdUsed[codim])
      {
        for(int i=0; i<en.template count<codim> (); i++)
        {
          lset.insert( hIndexSet. template subIndex<codim> (en,i) );   
        }
      }
    }
    
    static inline void removeCodims (const HSetImp & hIndexSet , 
        CodimLeafSet (&cls)[ncodim], const EntityType & en , bool (&cdUsed)[ncodim])
    {
      enum { codim = 1 };
      CodimLeafSet & lset = cls[codim];
      // if codim is already used, then also remove entities of this codim
      if(cdUsed[codim])
      {
        for(int i=0; i<en.template count<codim> (); i++)
        {
          lset.remove( hIndexSet. template subIndex<codim> (en,i) );   
        }
      }
    }
  };

  //! type of this class 
  typedef DefaultAdaptiveLeafIndexSet < GridType > ThisType;
  
  // my type, to be revised 
  enum { myType = 6 };

  typedef typename GridType :: Traits :: LeafIndexSet LeafIndexSetType;
  typedef typename GridType :: Traits :: LocalIdSet LocalIdSetType;
  typedef typename GridType :: template Codim<0> :: Entity EntityCodim0Type;
  const LeafIndexSetType & leafSet_; 
  const LocalIdSetType & idSet_; 

  enum { dim = GridType :: dimension };

  // flag for codim is in use or not 
  mutable bool codimUsed_ [ncodim];
  
  // true if all entities that we use are marked as USED 
  bool marked_;

  // true if the used entities were marked by grid walkthrough 
  bool markAllU_;

  // true if any of the higher codims is used 
  mutable bool higherCodims_;

  typedef PersistentLeafIndexSet<LocalIdSetType,LeafIndexSetType,0> PersistentLeafIndexSetType;
  mutable PersistentLeafIndexSetType * pLeafSet_;

public:
  //! type traits of this class
  typedef LeafIteratorTypes<GridType> Traits; 

  //! Constructor
  DefaultAdaptiveLeafIndexSet (const GridType & grid) 
    : DefaultGridIndexSetBase <GridType> (grid) ,  
    leafSet_( grid.leafIndexSet() ) , 
    idSet_  ( grid.localIdSet() ), 
    marked_ (false) , markAllU_ (false) , higherCodims_ (true), 
    pLeafSet_(0)
  {
    if(!pLeafSet_) pLeafSet_ = new PersistentLeafIndexSetType(idSet_,leafSet_); 
    // codim 0 is used by default
    codimUsed_[0] = true;

    // all higher codims are not used by default
    for(int i=1; i<ncodim; i++) codimUsed_[i] = false;
    //for(int i=1; i<ncodim; i++) codimUsed_[i] = true;

    // set the codim of each Codim Set. 
    //for(int i=0; i<ncodim; i++) codimLeafSet_[i].setCodim( i );
  }

  //! Destructor
  virtual ~DefaultAdaptiveLeafIndexSet () {};

  //! return type of index set, for GrapeDataIO
  int type () const { return myType; }

  //****************************************************************
  //
  //  INTERFACE METHODS for DUNE INDEX SETS 
  //
  //****************************************************************

  //! return global index 
  //! for dof mapper 
  // --index 
  template <class EntityType>
  int index (const EntityType & en) const
  {
    // this IndexWrapper provides specialisations for each codim 
    // see this class above 
    //return leafSet_.index(en); 
    return this->template index<0> (en,0);
  }
  
  //! return subIndex of given entity
  // see specialisation for codim 0 below 
  template <int cd>
  int subIndex (const EntityCodim0Type & en, int num) const
  {
    return this->template index<cd> (en,num);
    //return leafSet_.template subIndex<cd> (en,num);
  }

  //! return size of grid entities per level and codim 
  int size ( int codim , GeometryType type ) const
  {
    return persistentLeafSet().size();
    //return leafSet_.size(codim,type); 
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
  
  /** @brief Iterator to one past the last entity of given codim for partition type
   *  Here the grids leaf iterator is used 
   */
  template<int cd, PartitionIteratorType pitype>
  typename LeafIteratorTypes<GridType>::template Codim<cd>::
    template Partition<pitype>::Iterator end () const
  {
    return leafSet_.template end<cd,pitype> ();
  }

  /** @brief Iterator to first entity of given codimension and partition type.
   *  Here the grids leaf iterator is used 
   */
  template<int cd, PartitionIteratorType pitype>
  typename LeafIteratorTypes<GridType>::template Codim<cd>::
    template Partition<pitype>::Iterator begin () const
  {
    return leafSet_.template begin<cd,pitype> ();
  }

  //****************************************************************
  //
  //  METHODS for Adaptation with DofManger 
  //
  //****************************************************************

  //! insert index for father, mark childs index for removal  
  template <class EntityType>
  void restrictLocal ( EntityType &father, EntityType &son, bool initialize ) 
  {
    // important, first remove old, because 
    // on father indices might be used aswell 
    removeOldIndex( son );
    insertNewIndex( father );
  }

  //! insert indices for children , mark fathers index for removal  
  template <class EntityType>
  void prolongLocal ( EntityType &father, EntityType &son, bool initialize ) 
  {
    // important, first remove old, because 
    // on children indices might be used aswell 
    removeOldIndex( father );
    insertNewIndex( son );
  }

  //! insert new index to set 
  void insertNewIndex (const typename GridType::template Codim<0>::Entity & en )  
  {
    this->insert( en );
    marked_ = true;
  }

  //! Unregister entity which will be removed from the grid
  void removeOldIndex (const typename GridType::template Codim<0>::Entity & en )
  {
    this->remove( en ); 
    marked_ = true;
  }

  void resize () 
  {
    // give all entities that lie below the old entities new numbers 
    //markAllBelowOld ();
  }

  //! for dof manager, to check whether it has to copy dof or not 
  bool indexIsNew (int num, int codim) const
  {
    return persistentLeafSet().indexIsNew(num);
  }

  bool needsCompress () const { return true; }

  //! make to index numbers consecutive 
  //! return true, if at least one hole was closed 
  bool compress ()
  {
    //std::cout << "Start compressing index set\n";
    persistentLeafSet().resize();
    // if not marked, mark which indices are still used 
    markAllUsed(); 

    // true if a least one dof must be copied 
    bool haveToCopy = persistentLeafSet().compress(); 

    // next turn mark again 
    marked_   = false;
    markAllU_ = false;
    
    //std::cout << "Finished compressing index set\n";
    return haveToCopy;
  }

  //! memorise index 
  // --insert
  void insert (const EntityCodim0Type & en)
  {
    persistentLeafSet().insert(en);
  }

  //! set indices to unsed so that they are cleaned on compress  
  // --remove
  void remove (const EntityCodim0Type & en)
  {
    persistentLeafSet().remove(en);
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
  int index (const EntityType & en, int num) const
  {
    //return IndexWrapper<ThisType,LeafIndexSetType,EntityType,EntityType::codimension,codim>::
    //         index(*this,leafSet_,en,num);
    return persistentLeafSet().index(en,num);
  }
  
  //! return number of exisiting holes 
  int numberOfHoles(int codim)  const
  {
    return persistentLeafSet().numberOfHoles();
  }

  //! return old index, for dof manager only 
  int oldIndex (int num, int codim ) const
  {
    return persistentLeafSet().oldIndex(num);
  }

  //! return new index, for dof manager only returns index 
  int newIndex (int num , int codim ) const
  {
    return persistentLeafSet().newIndex(num);
  }

private:
  // insert index if entities lies below used entity, return 
  // false if not , otherwise return true
  bool insertNewIndex (const EntityCodim0Type & en, bool isLeaf , bool canInsert )
  {
    assert( pLeafSet_ );
    // if entity isLeaf then we insert index 
    if(isLeaf)
    {
      this->insert (en );
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
      //this->insert ( en );
      this->remove ( en );
      // set unused here, because index is only needed for prolongation 
    }
    return true;
  }

  PersistentLeafIndexSetType & persistentLeafSet() 
  { 
    assert( pLeafSet_ );
    return *pLeafSet_; 
  }
  const PersistentLeafIndexSetType & persistentLeafSet() const 
  { 
    assert( pLeafSet_ );
    return *pLeafSet_; 
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
    //std::cout << "Checked all old indices \n";
    marked_ = true;
  }
  
  //! give all entities that lie below the old entities new numbers 
  //! here we need the hierarchic iterator because for example for some
  //! grid more the one level of new elements can be created during adaption 
  //! there for we start to give new number for all elements below the old
  //! element 
  void markAllBelowOld () 
  {
    typedef typename GridType::template Codim<0>::LevelIterator LevelIteratorType; 

    int maxlevel = this->grid_.maxLevel();
    
    for(int level = 0; level<=maxlevel; level++)
    {
      LevelIteratorType levelend    = this->grid_.template lend  <0> (level);
      for(LevelIteratorType levelit = this->grid_.template lbegin<0> (level);
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
  template <int codim>
  int countElements( GeometryType type ) const 
  {
    typedef typename Traits :: template Codim <codim> :: 
        template Partition<All_Partition> :: Iterator IteratorType;

    int count = 0;
    IteratorType endit  = end<codim,All_Partition> ();
    for(IteratorType it = begin<codim,All_Partition> (); it != endit; ++it)
    {
      if( it->geometry().type() == type ) count++;
    }
    return count; 
  }
  
  // print interal data, for debugging only 
  // print if only done, if DEBUG_LEAFINDEXSET is defined 
  void print (const char * msg, bool oldtoo = false ) const;

public:

  // write indexset to xdr file 
  bool write_xdr(const std::basic_string<char> filename, int timestep) 
  {
    /*
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

    for(int i=0; i<ncodim; i++) codimLeafSet_[i].processXdr(&xdrs);

    xdr_destroy(&xdrs);
    fclose(file);
    */
    return true;
  }

  //! read index set from given xdr file 
  bool read_xdr(const std::basic_string<char> filename , int timestep)
  {
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
      for(int i=0; i<ncodim; i++) codimLeafSet_[i].processXdr(&xdrs);
    }

    xdr_destroy(&xdrs);
    fclose(file);
    */

    print("read Index set ");
    return true;
  }

}; // end of class AdaptiveLeafIndexSet 

template <class GridType> 
inline void DefaultAdaptiveLeafIndexSet<GridType>::
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


