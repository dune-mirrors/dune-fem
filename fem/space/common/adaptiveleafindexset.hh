#ifndef DUNE_LEAFINDEXSET_HH
#define DUNE_LEAFINDEXSET_HH

//- Dune includes 
#include <dune/grid/common/defaultindexsets.hh>

//- local includes 
#include <dune/fem/space/common/codimindexset.hh>

namespace Dune {

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
class AdaptiveLeafIndexSet : 
  public IndexSet<GridType, AdaptiveLeafIndexSet<GridType>, DefaultLeafIteratorTypes<GridType> >,
  public DefaultGridIndexSetBase <GridType>
{
public:
  enum { ncodim = GridType::dimension + 1 };

  enum INDEXSTATE { NEW, USED, UNUSED };
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
  template <class AdLeafSet, class HSetImp, class CodimLeafSet, class EntityType, int enCodim >
  struct DirectIndexWrapper
  {
    static inline int index (const AdLeafSet & ls , const HSetImp & hset, const CodimLeafSet &cls, 
                             const EntityType & en, bool cdUsed )
    {
      // if not setup for codim yet, do setup.
      if(!cdUsed) ls.template setUpCodimSet<enCodim> ();
      assert(cls.index ( hset.index( en ) ) >= 0 );
      return cls.index ( hset.index( en ) );
    }
  };

  // direct index return from method index (const EntityType & en )
  // for codim 0, do no checking
  template <class AdLeafSet, class HSetImp, class CodimLeafSet, class EntityType>
  struct DirectIndexWrapper<AdLeafSet,HSetImp,CodimLeafSet,EntityType,0>
  {
    static inline int index (const AdLeafSet & ls , const HSetImp & hset, const CodimLeafSet &cls, 
                             const EntityType & en ,  bool cdUsed )
    {
      assert(cls.index ( hset.index( en ) ) >= 0 );
      return cls.index ( hset.index( en ) );
    }
  };

  //**************************************************************//

  // index return from method index (const EntityType & en, int num )
  // this puts the index method and the subIndex methods together 
  template <class AdLeafSet, class HSetImp, class CodimLeafSet, class EntityType, int enCodim, int codim >
  struct IndexWrapper
  {
    static inline int index (const AdLeafSet & ls , const HSetImp & hset, const CodimLeafSet (&cls)[ncodim], 
                             const EntityType & en , int num , bool (&cdUsed)[ncodim] )
    {
      assert(cls[codim].index ( hset.index( en ) ) >= 0 );
      if(!cdUsed[codim]) ls.template setUpCodimSet<codim> ();
      return cls[codim].index ( hset.index( en ) );
    }
  };

  template <class AdLeafSet, class HSetImp, class CodimLeafSet, class EntityType>
  struct IndexWrapper<AdLeafSet,HSetImp,CodimLeafSet,EntityType,0,0>
  {
    static inline int index (const AdLeafSet & ls , const HSetImp & hset, const CodimLeafSet (&cls)[ncodim], 
                             const EntityType & en , int num ,  bool (&cdUsed)[ncodim] )
    {
      enum { codim = 0 };
      // check if we have index for given entity
      assert(cls[codim].index ( hset.index( en ) ) >= 0 );
      return cls[codim].index ( hset.index( en ) );
    }
  };

  //! if codim > codim of entity use subIndex 
  template <class AdLeafSet, class HSetImp, class CodimLeafSet, class EntityType>
  struct IndexWrapper<AdLeafSet,HSetImp,CodimLeafSet,EntityType,0,1>
  {
    static inline int index (const AdLeafSet & ls , const HSetImp & hset, const CodimLeafSet (&cls)[ncodim], 
                             const EntityType & en , int num ,  bool (&cdUsed)[ncodim] )
    {
      enum { codim = 1 };
      if(!cdUsed[codim]) ls.template setUpCodimSet<codim> ();
      assert(cls[codim].index ( hset.template subIndex<codim>( en , num ) ) >= 0 );
      return cls[codim].index ( hset.template subIndex<codim>( en , num ) );
    }
  };

  //! if codim > codim of entity use subIndex 
  template <class AdLeafSet, class HSetImp, class CodimLeafSet, class EntityType>
  struct IndexWrapper<AdLeafSet,HSetImp,CodimLeafSet,EntityType,0,2>
  {
    static inline int index (const AdLeafSet & ls , const HSetImp & hset, const CodimLeafSet (&cls)[ncodim], 
                             const EntityType & en , int num ,  bool (&cdUsed)[ncodim] )
    {
      enum { codim = 2 };
      assert( cls[codim].myCodim () == codim );
      
      if(!cdUsed[codim]) ls.template setUpCodimSet<codim> ();
      assert( cdUsed[codim] );
      assert(cls[codim].index ( hset.template subIndex<codim>( en , num ) ) >= 0 );
      return cls[codim].index ( hset.template subIndex<codim>( en , num ) );
    }
  };

  //! if codim > codim of entity use subIndex 
  template <class AdLeafSet, class HSetImp, class CodimLeafSet, class EntityType>
  struct IndexWrapper<AdLeafSet,HSetImp,CodimLeafSet,EntityType,0,3>
  {
    static inline int index (const AdLeafSet & ls , const HSetImp & hset, const CodimLeafSet (&cls)[ncodim], 
                             const EntityType & en , int num ,  bool (&cdUsed)[ncodim] )
    {
      enum { codim = 3 };
      if(!cdUsed[codim]) ls.template setUpCodimSet<codim> ();
      assert(cls[codim].index ( hset.template subIndex<codim>( en , num ) ) >= 0 );
      return cls[codim].index ( hset.template subIndex<codim>( en , num ) );
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
  typedef AdaptiveLeafIndexSet < GridType > ThisType;
  
  // my type, to be revised 
  enum { myType = 5 };

  typedef CodimIndexSet CodimIndexSetType; 
  mutable CodimIndexSetType codimLeafSet_[ncodim];

  // type of Hset Selector 
  typedef HierarchicIndexSetSelector<GridType> SelectorType;

  // my index set type 
  typedef typename SelectorType :: HierarchicIndexSet HIndexSetType;

  typedef typename GridType :: template Codim<0> :: Entity EntityCodim0Type;
  const HIndexSetType & hIndexSet_; 

  enum { dim = GridType :: dimension };

  // flag for codim is in use or not 
  mutable bool codimUsed_ [ncodim];
  
  // true if all entities that we use are marked as USED 
  bool marked_;

  // true if the used entities were marked by grid walkthrough 
  bool markAllU_;

  // true if any of the higher codims is used 
  mutable bool higherCodims_;

  mutable bool compressed_;

public:
  //! type traits of this class (see defaultindexsets.hh)
  typedef DefaultLeafIteratorTypes<GridType> Traits; 

  //! Constructor
  AdaptiveLeafIndexSet (const GridType & grid) 
    : DefaultGridIndexSetBase <GridType> (grid) ,  
    hIndexSet_( SelectorType::hierarchicIndexSet(grid) ) , 
    marked_ (false) , markAllU_ (false) , higherCodims_ (false) 
    //marked_ (false) , markAllU_ (false) , higherCodims_ (true) 
    , compressed_(true) // at start the set is compressed 
  {
    // codim 0 is used by default
    codimUsed_[0] = true;
    // all higher codims are not used by default
    for(int i=1; i<ncodim; i++) codimUsed_[i] = false;
    //for(int i=1; i<ncodim; i++) codimUsed_[i] = true;

    // set the codim of each Codim Set. 
    for(int i=0; i<ncodim; i++) codimLeafSet_[i].setCodim( i );

    resizeVectors();
    // give all entities that lie below the old entities new numbers 
    markAllUsed ();
  }

  //! Destructor
  virtual ~AdaptiveLeafIndexSet () {};

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
    enum { codim = EntityType::codimension };
    // this IndexWrapper provides specialisations for each codim 
    // see this class above 
    return DirectIndexWrapper<ThisType,HIndexSetType,CodimIndexSetType,EntityType,codim> ::
      index ( *this , hIndexSet_, codimLeafSet_[codim], en, codimUsed_[codim] );
  }
  
  //! return subIndex of given entity
  // see specialisation for codim 0 below 
  template <int cd>
  int subIndex (const EntityCodim0Type & en, int num) const
  {
    // this IndexWrapper provides specialisations for each codim 
    // see this class above 
    enum { enCodim = 0 }; // codim of entity is 0 
    return IndexWrapper<ThisType,HIndexSetType,CodimIndexSetType,EntityCodim0Type,enCodim, cd>::
             index(*this, hIndexSet_, codimLeafSet_, en , num , codimUsed_);
  }

  //! return size of grid entities per level and codim 
  int size (GeometryType type) const
  {
    int codim=GridType::dimension-type.dim();
    if( !codimUsed_[codim] )
    {
      assert( hIndexSet_.geomTypes(codim).size() == 1 ); 
      return CountElements<ThisType,dim>::count(*this,codim,type);
    }
    return codimLeafSet_[codim].size();
  }
  
  //! return size of grid entities per level and codim 
  int size ( int codim ) const
  {
    assert( hIndexSet_.geomTypes(codim).size() == 1 ); 
    return size(hIndexSet_.geomTypes(codim)[0]);
  }
  
  //! returns vector with geometry tpyes this index set has indices for
  const std::vector <GeometryType> & geomTypes (int codim) const 
  {
    return hIndexSet_.geomTypes(codim);
  }
  
  /** @brief Iterator to one past the last entity of given codim for partition type
   *  Here the grids leaf iterator is used 
   */
  template<int cd, PartitionIteratorType pitype>
  typename DefaultLeafIteratorTypes<GridType>::template Codim<cd>::
    template Partition<pitype>::Iterator end () const
  {
    return this->grid_.template leafend<cd,pitype> ();
  }

  /** @brief Iterator to first entity of given codimension and partition type.
   *  Here the grids leaf iterator is used 
   */
  template<int cd, PartitionIteratorType pitype>
  typename DefaultLeafIteratorTypes<GridType>::template Codim<cd>::
    template Partition<pitype>::Iterator begin () const
  {
    return this->grid_.template leafbegin<cd,pitype> ();
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
    // here we have to add the support of higher codims 
    resizeVectors();
    
    this->insert( en );
    marked_ = true;
  }

  //! Unregister entity which will be removed from the grid
  void removeOldIndex (const typename GridType::template Codim<0>::Entity & en )
  {
    this->remove( en ); 
    marked_ = true;
  }

  //! reallocate the vector for new size
  void resizeVectors()
  {
    codimLeafSet_[0].resize( hIndexSet_.size(0) );
    if(higherCodims_)
    {
      for(int i=1; i<ncodim; i++) 
      {  
        if(codimUsed_[i]) 
        {
          //std::cout << "resize codim " << i << "\n"; 
          codimLeafSet_[i].resize( hIndexSet_.size(i) );
        }
      }
    }
  }

  //! if grid has changed, resize index vectors, and create 
  //! indices for new entities, new entities are entities that 
  //! lie below the old entities 
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
      for(int i=1; i<ncodim; i++) 
        haveToCopy = (codimLeafSet_[i].compress()) ? true : haveToCopy;
    }

    // next turn mark again 
    marked_   = false;
    markAllU_ = false;
    
    compressed_ = true;
    return haveToCopy;
  }

  //! this index set needs compress after adaptation, here true is returned  
  bool needsCompress () const { return true; }

  //! memorise index 
  // --insert
  void insert (const EntityCodim0Type & en)
  {
    if( !codimLeafSet_[0].exsits( hIndexSet_.index(en) ) ) 
    {
      codimLeafSet_[0].insert ( hIndexSet_.index(en) );
      if(higherCodims_)
      {
        PartialSpec<HIndexSetType,CodimIndexSetType,EntityCodim0Type,dim> :: 
         iterateCodims ( hIndexSet_, codimLeafSet_, en , codimUsed_ ); 
      }
    }
    compressed_ = false;
  }

  //! set indices to unsed so that they are cleaned on compress  
  // --remove
  void remove (const EntityCodim0Type & en)
  {
    // if state is NEW or USED the index of all entities is removed 
    if( codimLeafSet_[0].exsits( hIndexSet_.index(en) ) ) 
    {
      codimLeafSet_[0].remove ( hIndexSet_.index(en) );
      if(higherCodims_)
      {
        PartialSpec<HIndexSetType,CodimIndexSetType,EntityCodim0Type,dim> :: 
          removeCodims ( hIndexSet_, codimLeafSet_, en , codimUsed_ ); 
      }
    }
    compressed_ = false;
  }

  //! return approximate size that is used during restriction 
  int additionalSizeEstimate () const 
  { 
    int addSize = 0; 
    for(int i=0; i<ncodim; i++) addSize += codimLeafSet_[i].additionalSizeEstimate(); 
    return addSize;
  }

  //! return global index 
  //! for dof mapper 
  // --index 
  template <int codim, class EntityType>
  int index (const EntityType & en, int num) const
  {
    return IndexWrapper<ThisType,HIndexSetType,CodimIndexSetType,EntityType,EntityType::codimension,codim>::
           index(*this,hIndexSet_,codimLeafSet_,en,num,codimUsed_);
  }
 
  //! return number of holes of the sets indices 
  int numberOfHoles ( int codim ) const
  {
    assert( codimUsed_[codim] );
    return codimLeafSet_[codim].numberOfHoles(); 
  }

  //! return old index, for dof manager only 
  int oldIndex (int num, int codim ) const
  {
    assert( codimUsed_[codim] );
    return codimLeafSet_[codim].oldIndex(num); 
  }

  //! return new index, for dof manager only returns index 
  int newIndex (int num , int codim ) const
  {
    assert( codimUsed_[codim] );
    return codimLeafSet_[codim].newIndex(num); 
  }

private:
  // insert index if entities lies below used entity, return 
  // false if not , otherwise return true
  bool insertNewIndex (const EntityCodim0Type & en, bool isLeaf , bool canInsert )
  {
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
      // if index >= 0, then all children may  also appear in the set 
      // from now on, indices can be inserted 
      if( codimLeafSet_[0].index( hIndexSet_.index (en ) ) >= 0 )
      {
        return true; 
      }

      // we have to go deeper 
      return false;
    }
    else 
    {
      // we insert to get an index 
      this->insert ( en );
      // we remove to unmark, because this is not a leaf entity 
      this->remove ( en );
    }
    return true;
  }

  //! mark indices that are still used and give new indices to 
  //! elements that need one 
  void markAllUsed () 
  {
    // unset all indices 
    for(int i=0; i<ncodim; i++) 
      if(codimUsed_[i]) codimLeafSet_[i].set2Unused(); 
    
    typedef typename GridType:: template Codim<0> :: LeafIterator LeafIteratorType; 
    // walk over leaf level on locate all needed entities  
    LeafIteratorType endit  = this->grid_.template leafend<0>   ();
    for(LeafIteratorType it = this->grid_.template leafbegin<0> (); 
        it != endit ; ++it )
    {
      this->insert( *it );
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
    
    typedef typename GridType:: template Codim<codim> :: LeafIterator LeafIteratorType; 
    // walk over leaf level on locate all needed entities  
    LeafIteratorType endit  = this->grid_.template leafend<codim>  ();
    for(LeafIteratorType it = this->grid_.template leafbegin<codim>(); 
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
    typedef typename GridType::template Codim<0>::LevelIterator LevelIteratorType; 

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
    return true;
  }

  //! read index set from given xdr file 
  bool read_xdr(const std::basic_string<char> filename , int timestep)
  {
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

    return true;
  }

}; // end of class AdaptiveLeafIndexSet 

} // end namespace Dune 
#endif
