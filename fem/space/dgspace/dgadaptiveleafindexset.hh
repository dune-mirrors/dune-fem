#ifndef DUNE_DGADAPTIVELEAFINDEXSET_HH
#define DUNE_DGADAPTIVELEAFINDEXSET_HH

//- Dune includes 
#include <dune/common/exceptions.hh>

//- local includes 
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/gridpart/dunefemindexsets.hh>
#include <dune/fem/gridpart/codimindexset.hh>

#include <dune/fem/io/file/xdrio.hh>

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
  public ConsecutivePersistentIndexSet<
        GridType, 
        DGAdaptiveLeafIndexSet<GridType>, 
        DefaultLeafIteratorTypes<GridType> 
          >
{
public:
  enum { ncodim = GridType::dimension + 1 };

protected:
  //! type of base class 
  typedef ConsecutivePersistentIndexSet<
        GridType, 
        DGAdaptiveLeafIndexSet<GridType>, 
        DefaultLeafIteratorTypes<GridType> 
          > BaseType;

public:  
  //! type of index 
  typedef typename BaseType :: IndexType IndexType;

protected:
  template <class CodimLeafSet,class HIndexSet, int codim>
  struct Index
  {
    template <class EntityType> 
    inline static int index(const CodimLeafSet & lset, 
                               const HIndexSet & hset, 
                               const EntityType & en) 
    {
      DUNE_THROW(NotImplemented,"DGAdaptiveLeafIndexSet does not support indices for higher codimension!");
      return 0;
    }
  };

  template <class CodimLeafSet,class HIndexSet>
  struct Index<CodimLeafSet,HIndexSet,0>
  {
    template <class EntityType> 
    inline static int index(const CodimLeafSet & lset, 
                               const HIndexSet & hset, 
                               const EntityType & en) 
    {
      return lset.index( hset.index( en ) );
    }
  };
  
  //! type of this class 
  typedef DGAdaptiveLeafIndexSet < GridType > ThisType;

  //! fro consecutive method 
  friend class Conversion< ThisType, EmptyIndexSet> ;
  
  //! is true if grid is structured grid 
  enum { StructuredGrid = ! Capabilities::IsUnstructured<GridType>::v };
  
  // my type, to be revised 
  enum { myType = (StructuredGrid) ? -1 : 665 };

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
  int sequence_;

  //! no copy allowed
  DGAdaptiveLeafIndexSet (const DGAdaptiveLeafIndexSet & ) ;
public:
  //! type traits of this class (see defaultindexsets.hh)
  typedef DefaultLeafIteratorTypes<GridType> Traits; 

  //! Constructor
  DGAdaptiveLeafIndexSet (const GridType & grid) 
    : BaseType(grid) 
    , codimLeafSet_( this->dofManager_.memoryFactor() )
    , hIndexSet_( SelectorType::hierarchicIndexSet(grid) ) 
    , compressed_(true) // at start the set is compressed 
    , sequence_(this->dofManager_.sequence())
  {
    // set the codim of this codim set, here always 0
    codimLeafSet_.setCodim( 0 );

    // setup all needed indices 
    setupIndexSet();
  }

  //! Destructor
  virtual ~DGAdaptiveLeafIndexSet () {}

  //! return type of index set, for GrapeDataIO
  static int type () { return myType; }

  //! return name of index set 
  std::string name () const { return "DGAdaptiveLeafIndexSet"; }

  //****************************************************************
  //
  //  INTERFACE METHODS for DUNE INDEX SETS 
  //
  //****************************************************************
  //! \brief return size of grid entities per level and codim 
  int size (GeometryType type) const
  {
    return codimLeafSet_.size();
  }
  
  //! \brief return size of grid entities per level and codim 
  int size ( int codim ) const
  {
    assert( hIndexSet_.geomTypes(0).size() == 1 ); 
    return size(hIndexSet_.geomTypes(0)[0]);
  }
  
  //! \brief returns vector with geometry tpyes this index set has indices for
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
  //! \brief insert new index to set 
  void insertEntity (const typename GridType::template Codim<0>::Entity & en )  
  {
    resizeVectors();
    
    // insert entity in set 
    this->insertIndex( en );
  }

  //! \brief un-register entity which will be removed from the grid
  void removeEntity (const typename GridType::template Codim<0>::Entity & en )
  {
    // only indices that are contained should be removed 
    assert( codimLeafSet_.index( hIndexSet_.index( en )) >= 0 );
    this->removeIndex( en ); 
  }

  //! \brief reallocate the vector for new size
  void resizeVectors()
  {
    codimLeafSet_.resize( hIndexSet_.size(0) );
  }

  /** \brief 
     if grid has changed, resize index vectors, and create 
     indices for new entities, new entities are entities that 
     lie below the old entities.
  */
  void resize () 
  {
    // adjust size of vectors 
    resizeVectors();

#if HAVE_MPI
    if( StructuredGrid  && 
        this->grid_.comm().size() > 1 )
    {
      codimLeafSet_.clear();

      // this should only be the case of YaspGrid
      markAllBelowOld<Interior_Partition>();
      markAllBelowOld<All_Partition>();
      compressed_ = true;
    }
    else 
#endif
    {
      // use a hierarchic walk to mark new elements 
      markAllBelowOld<All_Partition> ();

#if HAVE_MPI
      if( this->grid_.comm().size() > 1 )
      {
        // make sure that also ghosts have indices 
        markAllUsed<Ghost_Partition>();
      }
#endif
    }
  }

  //- --compress 
  /** \brief 
     make to index numbers consecutive 
     return true, if at least one hole existed  
  */
  bool compress ()
  {
    // reset list of holes in any case 
    codimLeafSet_.clearHoles(); 

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

    // mark all still needed indices 
    setupIndexSet();

    // true if a least one dof must be copied 
    bool haveToCopy = codimLeafSet_.compress(); 

    compressed_ = true;
    sequence_ = this->dofManager_.sequence();
    return haveToCopy;
  }

  /** \brief return global index for dof mapper */
  //- --index 
  template <int codim, class EntityType>
  inline IndexType indexImp (const EntityType & en, int num) const
  {
    return Index<CodimIndexSetType,HIndexSetType,codim>::index(codimLeafSet_,hIndexSet_,en);
  }
 
  //! \brief return number of holes of the sets indices 
  int numberOfHoles ( const int codim ) const
  {
    assert( codim == 0 );
    return codimLeafSet_.numberOfHoles(); 
  }

  //! \brief return old index, for dof manager only 
  int oldIndex (const int num, const int codim ) const
  {
    assert( codim == 0 );
    return codimLeafSet_.oldIndex(num); 
  }

  //! \brief return new index, for dof manager only returns index 
  int newIndex (const int num , const int codim ) const
  {
    assert( codim == 0 );
    return codimLeafSet_.newIndex(num); 
  }

protected:
  //! \brief memorise index 
  //- --insertIndex
  void insertIndex(const EntityCodim0Type & en)
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

  //! \brief set indices to unsed so that they are cleaned on compress  
  //- --removeIndex
  void removeIndex(const EntityCodim0Type & en)
  {
    const int idx = hIndexSet_.index(en);
    // if state is NEW or USED the index of all entities is removed 
    if( codimLeafSet_.exists( idx ) ) 
    {
      codimLeafSet_.remove ( idx );
      compressed_ = false;
    }
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
      // if index >= 0, then all children may  also appear in the set 
      // from now on, indices can be inserted, otherwise go deeper 
      return codimLeafSet_.validIndex( hIndexSet_.index (en ) );
    }
    else 
    {
      // we insert to get an index 
      this->insertIndex( en );
      // we remove to unmark, because this is not a leaf entity 
      this->removeIndex( en );
    }
    return true;
  }

  //! mark all indices of interest 
  void setupIndexSet ()
  {
    // for structured grids clear all information 
    // this in only done when setting up grids or after 
    // read of parallel data on serial grids 
    if( StructuredGrid )
    {
      // clear all information 
      codimLeafSet_.clear();
    }

#if HAVE_MPI
    // for YaspGrid we need all interior indices first 
    // so we can use SGrid for the visualization :(
    if( StructuredGrid  &&
        this->grid_.comm().size() > 1 )
    {
      // we should only get here for YaspGrid
      markAllUsed<Interior_Partition> ();
      markAllUsed<All_Partition>();
    }
    else
#endif
    {
      // give all entities that lie on the leaf level new numbers 
      markAllUsed<All_Partition> ();
    }
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
      this->insertIndex( *it );
    }
  }
  
  //! give all entities that lie below the old entities new numbers 
  //! here we need the hierarchic iterator because for example for some
  //! grid more the one level of new elements can be created during adaption 
  //! there for we start to give new number for all elements below the old
  //! element 
  template <PartitionIteratorType pitype>
  void markAllBelowOld() 
  {
    codimLeafSet_.set2Unused(); 

    typedef typename GridType::
      template Codim<0>::template Partition<pitype>:: LevelIterator LevelIteratorType;

    // iterate over macro level and check all entities hierachically 
    const LevelIteratorType macroend = this->grid_.template lend  <0,pitype> (0);
    for(LevelIteratorType macroit = this->grid_.template lbegin<0,pitype> (0);
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
    const bool isNew = insertNewIndex( en , en.isLeaf() , wasNew );

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

public:
  //! \brief write indexset to xdr file 
  bool write_xdr(const std::basic_string<char> filename, int timestep) 
  {
    const char *path = "";
    std::string fnstr = genFilename(path,filename, timestep);

    // create xdr write stream 
    XDRWriteStream xdr(fnstr);

    // check type of set 
    int type = myType;
    xdr.inout( type );

    return codimLeafSet_.processXdr(xdr);
  }

  //! \brief read index set from given xdr file 
  bool read_xdr(const std::string filename , int timestep)
  {
    const char *path = "";
    std::string fnstr = genFilename(path,filename, timestep);

    // create xdr read stream 
    XDRReadStream xdr(fnstr);

    // check type 
    int type = myType;
    xdr.inout( type );

    // check type 
    if( (type != 2) && (type != myType) )
    {
      DUNE_THROW(InvalidStateException,"DGAdaptiveLeafIndexSet::read_xdr: wrong type " << type << " given! Expected " << myType);
    }

    bool success = codimLeafSet_.processXdr(xdr);

    // in parallel runs we have to compress here
    if(this->grid_.comm().size() > 1) 
    {
      // mark for compress 
      compressed_ = false;
    }

    // for parallel data read in serial program we need compression 
    if(StructuredGrid && (hIndexSet_.size(0) != codimLeafSet_.size()) )
    {
      compressed_ = false;
    }
    return success;
  }
}; // end of class AdaptiveLeafIndexSet 

} // end namespace Dune 
#endif
