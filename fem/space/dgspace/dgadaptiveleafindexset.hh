#ifndef DUNE_DGADAPTIVELEAFINDEXSET_HH
#define DUNE_DGADAPTIVELEAFINDEXSET_HH

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

  Note that only codim = 0 is supported by this index set. 
*/
template <class GridType>
class DGAdaptiveLeafIndexSet : 
  public IndexSet<GridType, DGAdaptiveLeafIndexSet<GridType>, DefaultLeafIteratorTypes<GridType> >,
  public DefaultGridIndexSetBase <GridType>
{
public:
  enum { ncodim = GridType::dimension + 1 };

  template <class CodimLeafSet,class HIndexSet, int codim>
  struct SubIndex
  {
    template <class EntityType> 
    static int subIndex(CodimLeafSet & lset, const HIndexSet & hset, const EntityType & en) 
    {
      return 0;
    }
  };
  
  template <class CodimLeafSet,class HIndexSet>
  struct SubIndex<CodimLeafSet,HIndexSet,0>
  {
    template <class EntityType> 
    static int subIndex(CodimLeafSet & lset, const HIndexSet & hset, const EntityType & en) 
    {
      return lset.index( hset.index( en ) );
    }
  };
  
private:
  //! type of this class 
  typedef DGAdaptiveLeafIndexSet < GridType > ThisType;
  
  // my type, to be revised 
  enum { myType = 665 };

  typedef CodimIndexSet CodimIndexSetType;
  mutable CodimIndexSetType codimLeafSet_;

  // type of Hset Selector 
  typedef HierarchicIndexSetSelector<GridType> SelectorType;

  // my index set type 
  typedef typename SelectorType :: HierarchicIndexSet HIndexSetType;

  typedef typename GridType :: template Codim<0> :: Entity EntityCodim0Type;
  const HIndexSetType & hIndexSet_; 

  enum { dim = GridType :: dimension };

  //! true if set is consecutive without any holes 
  mutable bool compressed_;

  typedef DofManager<GridType> DofManagerType; 
  typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;
  
public:
  static ThisType & instance (const GridType & grid) 
  {
    static ThisType set(grid);
    return set;
  }
  
  //! type traits of this class (see defaultindexsets.hh)
  typedef DefaultLeafIteratorTypes<GridType> Traits; 

  //! Constructor
  DGAdaptiveLeafIndexSet (const GridType & grid) 
    : DefaultGridIndexSetBase <GridType> (grid) 
    , codimLeafSet_( DofManagerFactoryType::getDofManager(grid).memoryFactor() )
    , hIndexSet_( SelectorType::hierarchicIndexSet(grid) ) 
    , compressed_(true) // at start the set is compressed 
  {
    // set the codim of this codim set, here always 0
    codimLeafSet_.setCodim( 0 );

    // give all entities that lie below the old entities new numbers 
    markAllUsed<All_Partition> ();
  }

  //! Destructor
  virtual ~DGAdaptiveLeafIndexSet () {};

  //! return type of index set, for GrapeDataIO
  static int type () { return myType; }

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
    assert( codim == 0 ); // only set for codim 0 
    return codimLeafSet_.index( hIndexSet_.index( en ) ) ;
  }
  
  //! return subIndex of given entity
  // see specialisation for codim 0 below 
  template <int cd>
  int subIndex (const EntityCodim0Type & en, int) const
  {
    // this IndexWrapper provides specialisations for each codim 
    // see this class above 
    return SubIndex<CodimIndexSetType,HIndexSetType,cd>::subIndex(codimLeafSet_,hIndexSet_,en);
  }

  //! return size of grid entities per level and codim 
  int size (GeometryType type) const
  {
    return codimLeafSet_.size();
  }
  
  //! return size of grid entities per level and codim 
  int size ( int codim ) const
  {
    assert( hIndexSet_.geomTypes(0).size() == 1 ); 
    return size(hIndexSet_.geomTypes(0)[0]);
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

  //! \brief returns true if entity is contained in index set 
  template <class EntityType>
  bool contains(const EntityType & en) const
  {
    enum { codim = EntityType::codimension };
    return (codim == 0) ? codimLeafSet_.exists( hIndexSet_.index( en ) ) : false;
  }
  
  //****************************************************************
  //
  //  METHODS for Adaptation with DofManger 
  //
  //****************************************************************

  //! insert index for father, mark childs index for removal  
  template <class EntityType>
  void restrictLocal ( EntityType& father, EntityType& son, bool initialize ) const
  {
    // important, first remove old, because 
    // on father indices might be used aswell 
    removeOldIndex( son );
    insertNewIndex( father );
  }

  //! insert indices for children , mark fathers index for removal  
  template <class EntityType>
  void prolongLocal ( EntityType& father, EntityType& son, bool initialize ) const
  {
    // important, first remove old, because 
    // on children indices might be used aswell 
    removeOldIndex( father );
    insertNewIndex( son );
  }
 
  //! insert new index to set 
  void insertNewIndex (const typename GridType::template Codim<0>::Entity & en )  
  {
    resizeVectors();
    
    // insert entity in set 
    this->insert( en );
  }

  //! Unregister entity which will be removed from the grid
  void removeOldIndex (const typename GridType::template Codim<0>::Entity & en )
  {
    // only indices that are contained should be removed 
    assert( codimLeafSet_.index( hIndexSet_.index( en )) >= 0 );
    this->remove( en ); 
  }

  //! reallocate the vector for new size
  void resizeVectors()
  {
    codimLeafSet_.resize( hIndexSet_.size(0) );
  }

  //! if grid has changed, resize index vectors, and create 
  //! indices for new entities, new entities are entities that 
  //! lie below the old entities 
  void resize () 
  {
    // adjust size of vectors 
    resizeVectors();

    // mark all children 
    markAllBelowOld<All_Partition> ();
  }

  //! this index set can be used for adaptive calculations 
  bool adaptive () const { return true; }

  //! make to index numbers consecutive 
  //! return true, if at least one hole existed  
  //- --compress 
  bool compress ()
  {
    // if set already compress, do noting 
    if(compressed_) return false;

    // mark all leaf elements  
    // needs a leaf traversal 
    markAllUsed<All_Partition> (); 

    // true if a least one dof must be copied 
    bool haveToCopy = codimLeafSet_.compress(); 

    compressed_ = true;
    return haveToCopy;
  }

  //! this index set needs compress after adaptation, here true is returned  
  bool needsCompress () const { return true; }

  //! memorise index 
  //- --insert
  void insert (const EntityCodim0Type & en)
  {
    const int idx = hIndexSet_.index(en);
    if( !codimLeafSet_.exists( idx ) ) 
    {
#if HAVE_MPI 
      // we need special treatment for ghosts 
      // ghosts should not be inlcuded in holes list 
      if(en.partitionType() == GhostEntity) 
      {
        codimLeafSet_.insertGhost ( idx );
      }
      else   
#endif
      {
        codimLeafSet_.insert ( idx );
      }
      compressed_ = false;
    }

    assert( codimLeafSet_.exists( idx ) );
  }

  //! set indices to unsed so that they are cleaned on compress  
  //- --remove
  void remove (const EntityCodim0Type & en)
  {
    const int idx = hIndexSet_.index(en);
    // if state is NEW or USED the index of all entities is removed 
    if( codimLeafSet_.exists( idx ) ) 
    {
      codimLeafSet_.remove ( idx );
      compressed_ = false;
    }
  }

  //! return approximate size that is used during restriction 
  int additionalSizeEstimate () const 
  {
    return codimLeafSet_.additionalSizeEstimate();
  }

  //! return global index 
  //! for dof mapper 
  //- --index 
  template <int codim, class EntityType>
  int index (const EntityType & en, int num) const
  {
    assert( codim == 0 );
    // make sure that requested index was generated before 
    assert( codimLeafSet_.index( hIndexSet_.index( en )) >= 0 );
    return codimLeafSet_.index( hIndexSet_.index( en )); 
  }
 
  //! return number of holes of the sets indices 
  int numberOfHoles ( int codim ) const
  {
    assert( codim == 0 );
    return codimLeafSet_.numberOfHoles(); 
  }

  //! return old index, for dof manager only 
  int oldIndex (int num, int codim ) const
  {
    assert( codim == 0 );
    return codimLeafSet_.oldIndex(num); 
  }

  //! return new index, for dof manager only returns index 
  int newIndex (int num , int codim ) const
  {
    assert( codim == 0 );
    return codimLeafSet_.newIndex(num); 
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
      if( codimLeafSet_.index( hIndexSet_.index (en ) ) >= 0 )
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
  template <PartitionIteratorType pitype>
  void markAllUsed () 
  {
    // make correct size of vectors 
    resizeVectors();
    
    // unset all indices 
    codimLeafSet_.set2Unused(); 
    
    typedef typename GridType:: template Codim<0> :: template 
        Partition<pitype> :: LeafIterator LeafIteratorType; 
    // walk over leaf level and locate all needed entities  
    LeafIteratorType endit  = this->grid_.template leafend<0,pitype>   ();
    for(LeafIteratorType it = this->grid_.template leafbegin<0,pitype> (); 
        it != endit ; ++it )
    {
      this->insert( *it );
    }
  }
  
  //! give all entities that lie below the old entities new numbers 
  //! here we need the hierarchic iterator because for example for some
  //! grid more the one level of new elements can be created during adaption 
  //! there for we start to give new number for all elements below the old
  //! element 
  template <PartitionIteratorType pitype>
  void markAllBelowOld () 
  {
    typedef typename GridType::template Codim<0>::
          template Partition<pitype> :: LevelIterator LevelIteratorType; 

    int maxlevel = this->grid_.maxLevel();
   
    codimLeafSet_.set2Unused(); 
    
    for(int level = 0; level<=maxlevel; ++level)
    {
      LevelIteratorType levelend    = this->grid_.template lend  <0,pitype> (level);
      for(LevelIteratorType levelit = this->grid_.template lbegin<0,pitype> (level);
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
  }

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

    codimLeafSet_.processXdr(&xdrs);

    xdr_destroy(&xdrs);
    fclose(file);
    return true;
  }

  //! read index set from given xdr file 
  bool read_xdr(const std::string filename , int timestep)
  {
    FILE   *file;
    XDR     xdrs;
    const char *path = "";

    std::string fnstr = genFilename(path,filename, timestep);
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

    codimLeafSet_.processXdr(&xdrs);

    xdr_destroy(&xdrs);
    fclose(file);

    // in parallel runs we have to compress here
    if(this->grid_.comm().size() > 1) 
    {
      compressed_ = false;
      compress();
    }

    return true;
  }

}; // end of class AdaptiveLeafIndexSet 

} // end namespace Dune 
#endif
