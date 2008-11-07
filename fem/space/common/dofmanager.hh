#ifndef DUNE_DOFMANAGER_HH
#define DUNE_DOFMANAGER_HH

//- System includes 
#include <cassert>
#include <string>
#include <list>

//- Dune includes 
#include <dune/common/stdstreams.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/interfaces.hh>

// here are the default grid index set defined 
#include <dune/fem/space/common/restrictprolonginterface.hh>
#include <dune/fem/storage/singletonlist.hh>

#include <dune/fem/io/parameter.hh>

//- local includes 
#include "dofmapper.hh"
#include "datacollector.hh"
#include "arrays.hh"

namespace Dune
{

/** @addtogroup DofManager  

    @{
**/

// forward declaration 
template <class GridType> class DofManager;
template <class DofManagerImp> class DofManagerFactory;

/** \brief SpecialArrayFeatures is a wrapper class to extend some array
    classes with some special features needed for the MemObject.
    There exsist a specialization for MutableArray. 
 */ 
template<class ArrayType>
struct SpecialArrayFeatures
{
  /** \brief value type of array, i.e. double */ 
  typedef typename ArrayType :: value_type ValueType;

  /** \brief return used memory size of Array */
  static size_t used(const ArrayType & array)  
  {
    return array.size() * sizeof(ValueType);
  }
  
  /** \brief set memory overestimate factor, here does nothing */ 
  static void setMemoryFactor(ArrayType & array, const double memFactor) 
  {
  }

  /** \brief move memory blocks backwards */
  static void memMoveBackward(ArrayType& array, const int length,
            const int oldStartIdx, const int newStartIdx)
  {
    // get new end of block which is offSet + (length of block - 1) 
    int newIdx = newStartIdx + length - 1; 
    // copy all entries backwards 
    for(int oldIdx = oldStartIdx+length-1; oldIdx >= oldStartIdx; --oldIdx, --newIdx )
    {
      // move value to new location 
      array[newIdx] = array[oldIdx];
    }
  }

  /** \brief move memory blocks forward */
  static void memMoveForward(ArrayType& array, const int length,
            const int oldStartIdx, const int newStartIdx)
  {
    const int upperBound = oldStartIdx + length;
    // get new off set that should be smaller then old one
    int newIdx = newStartIdx;
    for(int oldIdx = oldStartIdx; oldIdx<upperBound; ++oldIdx, ++newIdx )
    {
      // copy to new location 
      array[newIdx] = array[oldIdx];
    }
  }
};

///////////////////////////////////////////////////////////////////////
//
//  ManagedIndexSetInterface 
//
///////////////////////////////////////////////////////////////////////
/*! This class is the virtual interface for the index sets managed by
  the DofManager. The derived classes are of the type ManagedIndexSet<IndexSet>.
  This means we don't have to inherit every index set we want to use with
  this DofManager. 
  */
class ManagedIndexSetInterface
{
  ManagedIndexSetInterface(const ManagedIndexSetInterface& org);
protected:
  // use address of object as id 
  typedef const void * IdentifierType;
  // pointer to compare index sets 
  IdentifierType setPtr_;
  // reference counter 
  size_t referenceCounter_;
  
  template <class IndexSetType>
  ManagedIndexSetInterface(const IndexSetType& set) 
    : setPtr_( getIdentifier(set) ) , referenceCounter_(1) 
  {}

public:
  virtual ~ManagedIndexSetInterface () {}
  //! resize of index set 
  virtual void resize () = 0; 
  //! compress of index set 
  virtual bool compress () = 0;
  //! read and write method of index sets 
  virtual void read_xdr(const char * filename, int timestep) = 0;
  //! read and write method of index sets 
  virtual void write_xdr(const char * filename, int timestep) const = 0;

  //! increase reference counter 
  template <class IndexSetType>
  bool increaseReference(const IndexSetType& set) 
  {
    // if index sets are the same, return true and increase counter 
    return ( equals( set ) ) ? (++referenceCounter_,true) : false;
  } 

  template <class IndexSetType>
  bool decreaseReference(const IndexSetType& set) 
  {
    // if index sets are the same, 
    // decrease and return true if counter zero 
    return ( equals( set ) ) ? (--referenceCounter_ == 0) : false;
  }
  
private:
  template <class IndexSetType>
  bool equals(const IndexSetType& set) const
  {
    // if index sets are the same 
    return (getIdentifier(set) == setPtr_ );  
  }

  template <class IndexSetType>
  IdentifierType getIdentifier(const IndexSetType& set) const 
  {
    return static_cast<IdentifierType> (&set); 
  } 
    
};

template <class IndexSetType, class EntityType> class RemoveIndicesFromSet;
template <class IndexSetType, class EntityType> class InsertIndicesToSet;

template <class IndexSetType, class EntityType>
class ManagedIndexSet : public ManagedIndexSetInterface ,
        public LocalInlinePlus < ManagedIndexSet<IndexSetType,EntityType> , EntityType >
{
  typedef LocalInterface<EntityType> LocalIndexSetObjectsType;
protected: 
  // the dof set stores number of dofs on entity for each codim
  IndexSetType & indexSet_;

  // insertion and removal of indices 
  InsertIndicesToSet   <IndexSetType,EntityType> insertIdxObj_;
  RemoveIndicesFromSet <IndexSetType,EntityType> removeIdxObj_;

  LocalIndexSetObjectsType & indexSetList_; 
  LocalIndexSetObjectsType & insertList_; 
  LocalIndexSetObjectsType & removeList_; 

public:  
  //! type of base class 
  typedef ManagedIndexSetInterface BaseType;
  
  //! Constructor of MemObject, only to call from DofManager 
  ManagedIndexSet ( const IndexSetType & iset 
      , LocalIndexSetObjectsType & indexSetList
      , LocalIndexSetObjectsType & insertList 
      , LocalIndexSetObjectsType & removeList) 
   : BaseType( iset )
   , indexSet_ (const_cast<IndexSetType &> (iset)) 
   , insertIdxObj_(indexSet_), removeIdxObj_(indexSet_) 
   , indexSetList_(indexSetList) 
   , insertList_(insertList) 
   , removeList_(removeList)
  {
    this->setPtr_ = (void *) &indexSet_;
    
    indexSetList_ += *this;
    if( indexSet_.consecutive() ) 
    {
      insertList_ += insertIdxObj_; 
      removeList_ += removeIdxObj_;
    }
  } 

  //! desctructor 
  ~ManagedIndexSet () 
  {
    indexSetList_.remove( *this );
    if( indexSet_.consecutive() ) 
    {
      insertList_.remove( insertIdxObj_ ); 
      removeList_.remove( removeIdxObj_ );
    }
  }

  //! wrap resize of index set 
  void resize () 
  {
    indexSet_.resize();
  }
  
  //! wrap compress of index set 
  bool compress () 
  { 
    return indexSet_.compress(); 
  }

  //! call read_xdr of index set 
  virtual void read_xdr(const char * filename, int timestep)
  {
    indexSet_.read_xdr(filename,timestep); 
  }
  
  //! call write_xdr of index set 
  virtual void write_xdr(const char * filename, int timestep) const
  {
    indexSet_.write_xdr(filename,timestep);
  }
};

/////////////////////////////////////////////////////////////
//
// DofStorageInterface 
//
/////////////////////////////////////////////////////////////
/** \brief Interface class for a dof storage object to be stored in
    discrete functions */
class DofStorageInterface
{
protected:
  //! do not allow to create explicit instances 
  DofStorageInterface() {}

public:
  //! destructor 
  virtual ~DofStorageInterface() {};

  //! enable dof compression for dof storage (default is empty)
  virtual void enableDofCompression() { };

  //! returns name of dof storage 
  virtual const std::string& name () const  = 0;

  //! size of space, i.e. mapper.size()
  virtual int size () const = 0;
};


/** \brief Interface class for a dof storage object that can be managed
    (resized and compressed) by the DofManager 
*/
class ManagedDofStorageInterface : public DofStorageInterface
{
protected:
  //! do not allow to create explicit instances 
  ManagedDofStorageInterface() {}

public:
  //! destructor 
  virtual ~ManagedDofStorageInterface() {};

  //! resize memory 
  virtual void resize () = 0;
  //! resize memory 
  virtual void resize (int newSize) = 0;
  //! resize memory 
  virtual void reserve (int newSize) = 0;
  //! new size of space, i.e. after grid adaptation
  virtual int newSize () const = 0;
  //! compressed the underlying dof vector 
  virtual void dofCompress () = 0;
  //! return size of mem used by MemObject 
  virtual size_t usedMemorySize() const = 0;
};


template <class MemObjectType> class ResizeMemoryObjects;
template <class MemObjectType> class ReserveMemoryObjects;

/*! 
  A ManagedDofStorage holds the memory for one DiscreteFunction and the
  corresponding DofArrayMemory. If a DiscreteFunction is signed in by a
  function space, then such a MemObject is created by the DofManager. 
  The MemObject also knows the DofMapper from the function space which the
  discrete function belongs to. Here we dont know the exact type of the dof
  mapper therefore the methods newSize and calcInsertPoints of the mappers
  have to be virtual. This isnt a problem because this methods should only
  be called during memory reorganizing which is only once per timestep. 
*/
template <class GridImp, class MapperType , class DofArrayType>
class ManagedDofStorage : public ManagedDofStorageInterface
{
  // interface for MemObject lists
  typedef LocalInterface< int > MemObjectCheckType;
private:
  // type of this class 
  typedef ManagedDofStorage <GridImp, MapperType , DofArrayType> ThisType;

  typedef DofManager<GridImp> DofManagerType;
  typedef DofManagerFactory< DofManagerType > DofManagerFactoryType;

  // reference to dof manager 
  DofManagerType& dm_;

  // the dof set stores number of dofs on entity for each codim
  mutable MapperType &mapper_;

  // Array which the dofs are stored in 
  DofArrayType array_;

  // name of mem object, i.e. name of discrete function 
  std::string name_;

  typedef ResizeMemoryObjects < ThisType > ResizeMemoryObjectType;
  typedef ReserveMemoryObjects  < ThisType > ReserveMemoryObjectType;
  ResizeMemoryObjectType  resizeMemObj_;
  ReserveMemoryObjectType reserveMemObj_;

  // true if data need to be compressed 
  bool dataNeedCompress_;

  // prohibit copying 
  ManagedDofStorage(const ManagedDofStorage& );
public:
  //! Constructor of MemObject, only to call from DofManager 
  ManagedDofStorage ( const GridImp& grid,
                      const MapperType& mapper,
                      const std::string& name
                    )
    : dm_( DofManagerFactoryType :: getDofManager( grid ) ),
      mapper_ ( const_cast<MapperType& >(mapper)),
      array_( mapper_.size() ),
      name_ (name),
      resizeMemObj_(*this),
      reserveMemObj_(*this),
      dataNeedCompress_(false)
  {
    // add to dof manager 
    dm_.addDofStorage( *this );

    // set memory over estimate factor, only for DofArray 
    SpecialArrayFeatures<DofArrayType>::setMemoryFactor(array_,dm_.memoryFactor());
  }

  //! \brief destructor deleting MemObject from resize and reserve List
  ~ManagedDofStorage()
  {
    // remove from dof manager 
    dm_.removeDofStorage( *this );
  }

  //! return object that calls resize of this memory object 
  ResizeMemoryObjectType& resizeMemoryObject() { return resizeMemObj_; }

  //! return object that calls reserve of this memory object 
  ReserveMemoryObjectType& reserveMemoryObject() { return reserveMemObj_; }

  //! returns name of this vector 
  const std::string& name () const { return name_; }

  //! if grid changed, then calulate new size of dofset 
  inline int newSize () const
  {
    return mapper().newSize();
  }

  //! return size of underlying array 
  int size () const { return array_.size(); }

  //! resize the memory with the new size 
  void resize () 
  {
    resize ( newSize() );
  }

  //! resize the memory with the new size 
  void resize ( const int nSize ) 
  {
    // if the size is allready correct, do nothing
    if( nSize == size() ) return ;

    // resize memory 
    array_.resize( nSize );
    // update mapper and move array 
    moveToRear();
  }

  //! reserve memory for what is comming 
  inline void reserve ( const int needed )
  {
    // if index set is compressible, then add requested size 
    if( mapper().consecutive() )
    {
      const int nSize = size() + (needed * mapper().maxNumDofs());
      array_.reserve( nSize );
    }
    else 
    {
      // if compress is not needed just resize with given size 
      // therefore use newSize to enleage array 
      assert( ! mapper().consecutive() );
      array_.resize( newSize() );
    }
  }

  //! copy the dof from the rear section of the vector to the holes 
  void dofCompress () 
  {
    const int nSize = newSize();
    if( dataNeedCompress_ && mapper().consecutive() )
    {
      const int oldSize = mapper().size();
      // update mapper to new sizes 
      mapper().update();

      const int numBlocks = mapper().numBlocks();
      for( int block = 0; block < numBlocks; ++block )
      {
        moveToFront( oldSize, block );

        // run over all holes and copy array vules to new place 
        const int holes = mapper().numberOfHoles( block );
        for( int i = 0; i < holes; ++i )
        {
          const int oldIndex = mapper().oldIndex( i, block );
          const int newIndex = mapper().newIndex( i, block );

          assert( newIndex < nSize );
          array_[ newIndex ] = array_[ oldIndex ];
        }
      }
    }

    // store new size, which should be smaller then actual size 
    array_.resize( nSize );
  }
 
  //! return used memory size 
  size_t usedMemorySize() const 
  {
    return ((size_t) sizeof(ThisType) + SpecialArrayFeatures<DofArrayType>::used(array_)); 
  }

  //! enable dof compression for this MemObject
  void enableDofCompression() 
  {
    dataNeedCompress_ = true;
  }

  //! return reference to array for DiscreteFunction 
  DofArrayType & getArray() { return array_; } 

protected:
  inline MapperType &mapper () const
  {
    return mapper_;
  }
  
  // move array to rear insertion points 
  void moveToRear ()
  {
    // calculate new insertion points 
    const int oldSize = mapper().size();

    // update mapper to new sizes 
    mapper().update();

    // now check all blocks beginning with the largest 
    const int numBlocks = mapper().numBlocks();
    for( int block = numBlocks-1; block >= 0; --block )
    {
      // get old off set 
      const int oldOffSet = mapper().oldOffSet( block );
      // if off set is not zero  
      if( oldOffSet > 0 )
      {
        const int newOffSet = mapper().offSet( block );
        // get upperBound
        const int upperBound
          = (block == numBlocks - 1) ? oldSize : mapper().oldOffSet( block + 1 );
        const int blockSize = upperBound - oldOffSet;

        // move block backward 
        SpecialArrayFeatures< DofArrayType >
          :: memMoveBackward( array_, blockSize, oldOffSet, newOffSet );
      }
    }
  }

  //! move block to front again 
  void moveToFront ( const int oldSize, const int block )
  {
    // get insertion point from block
    const int oldOffSet = mapper().oldOffSet( block );

    const int numBlocks = mapper().numBlocks();
    // only if block is not starting from zero 
    if( oldOffSet > 0 )
    {
      const int newOffSet = mapper().offSet( block );
      
      // for last section upperBound is size 
      const int upperBound
        = (block == numBlocks - 1) ? oldSize : mapper().oldOffSet( block + 1 );
      const int blockSize = upperBound - oldOffSet;

      // move block forward 
      SpecialArrayFeatures< DofArrayType >
        :: memMoveForward( array_, blockSize, oldOffSet, newOffSet ); 
    }
  }
};

//! default implementation for creating a managed dof storage 
template <class GridType, class MapperType, class DofStorageType>
static inline std::pair< DofStorageInterface* , DofStorageType* >
  allocateManagedDofStorage( const GridType& grid,
                             const MapperType& mapper,
                             const std::string& name,
                             const DofStorageType* )
{
  // create managed dof storage 
  typedef ManagedDofStorage< GridType, MapperType,
                             DofStorageType > ManagedDofStorageType;
      
  ManagedDofStorageType* mds =
    new ManagedDofStorageType( grid, mapper, name );
  assert( mds );

  // return pair with dof storage pointer and array pointer 
  return std::pair< DofStorageInterface* , DofStorageType* >
          ( mds , & mds->getArray () );
}



///////////////////////////////////////////////////////////////
//
//  RestrictPorlong for Index Sets 
//
///////////////////////////////////////////////////////////////

template <class IndexSetType, class EntityType>
class RemoveIndicesFromSet 
: public LocalInlinePlus < RemoveIndicesFromSet<IndexSetType,EntityType> , EntityType >
{
private:
  // the dof set stores number of dofs on entity for each codim
  IndexSetType & indexSet_;

public:  
  // Constructor of MemObject, only to call from DofManager 
  RemoveIndicesFromSet ( IndexSetType & iset ) : indexSet_ (iset) {} 

  //! apply wraps the removeEntity Method of the index set 
  inline void apply ( EntityType & en )
  {
    indexSet_.removeEntity( en );
  }
};

template <class IndexSetType, class EntityType>
class InsertIndicesToSet 
: public LocalInlinePlus < InsertIndicesToSet<IndexSetType,EntityType> , EntityType >
{
private:
  // the dof set stores number of dofs on entity for each codim
  IndexSetType & indexSet_;

public:  
  // Constructor of MemObject, only to call from DofManager 
  InsertIndicesToSet ( IndexSetType & iset ) : indexSet_ (iset) {} 

  //! apply wraps the insertEntity method of the index set
  inline void apply ( EntityType & en )
  {
    indexSet_.insertEntity( en );
  }
};

template <class MemObjectType> 
class ResizeMemoryObjects 
: public LocalInlinePlus < ResizeMemoryObjects < MemObjectType > , int > 
{
private:
  // the dof set stores number of dofs on entity for each codim
  MemObjectType & memobj_;

public:  
  // Constructor of MemObject, only to call from DofManager 
  ResizeMemoryObjects ( MemObjectType & mo ) : memobj_ (mo) {} 
  ResizeMemoryObjects ( const ResizeMemoryObjects& org ) 
    : memobj_(org.memobj_)
  {}

  // resize mem object, parameter not needed 
  inline void apply ( int & )
  {
    memobj_.resize();
  }
};

// this class is the object for a single MemObject to 
template <class MemObjectType> 
class ReserveMemoryObjects  
: public LocalInlinePlus < ReserveMemoryObjects < MemObjectType > , int > 
{
private:
  // the dof set stores number of dofs on entity for each codim
  MemObjectType & memobj_;

public:  
  // Constructor of MemObject, only to call from DofManager 
  ReserveMemoryObjects ( MemObjectType & mo ) : memobj_ (mo) {} 

  // reserve for at least chunkSize new values 
  inline void apply ( int & chunkSize )
  {
    memobj_.reserve( chunkSize );  
  }
};

// this is the dofmanagers object which is being used during restriction 
// and prolongation process for adding and removing indices to and from
// index sets which belong to functions that belong to that dofmanager
template <class DofManagerType , class RestrictProlongIndexSetType> 
class IndexSetRestrictProlong  : 
  public RestrictProlongInterface<RestrictProlongTraits<IndexSetRestrictProlong<DofManagerType,RestrictProlongIndexSetType> > >
{
  DofManagerType & dm_;
  
  RestrictProlongIndexSetType & insert_;
  RestrictProlongIndexSetType & remove_;
public: 
  
  IndexSetRestrictProlong ( DofManagerType & dm , RestrictProlongIndexSetType & is, RestrictProlongIndexSetType & rm ) 
    : dm_(dm) , insert_(is), remove_(rm) {}

  // required for interface
  typedef double RangeFieldType;
  void setFatherChildWeight (const RangeFieldType& val) const {
  }

  //! restrict data to father 
  template <class EntityType>
  inline void restrictLocal ( EntityType & father, EntityType & son , bool initialize ) const
  {
    insert_.apply( father );
    remove_.apply( son );

    // resize memory 
    dm_.resizeMemory();
  }  

  //! prolong data to children 
  template <class EntityType>
  inline void prolongLocal ( EntityType & father, EntityType & son , bool initialize ) const
  {
    remove_.apply( father );
    insert_.apply( son );
    
    // resize memory 
    dm_.resizeMemory();
  }
  
};

class DofManError : public Exception {};

/*! 
 The DofManager is responsible for managing memory allocation and freeing
 for all discrete functions living on the grid the manager belongs to. 
 There is only one DofManager per grid.
 Each discrete function knows its dofmanager and can sign in. 
 If the grid is adapted, then the
 dofmanager reorganizes the memory if necessary. The DofManager holds a
 list of MemObjects which manage the memory and the corresponding
 mapper so they can determine the size of new memory. 
 Furthermore the DofManager holds an IndexSet which the DofMapper needs for
 calculating the indices in the dof vector for a given entity and local dof
 number. This IndexSet is delivered to the mapper when a function space is
 created. The default value for the IndexSet is the DefaultIndexSet class
 which is mostly a wrapper for the grid indices. 
*/
// --DofManager 
template< class Grid > 
class DofManager
: public IsDofManager 
{
  typedef DofManager< Grid > ThisType;

  friend class DefaultSingletonFactory< const Grid*, ThisType >;
  friend class DofManagerFactory< ThisType >;
  friend class Conversion< ThisType, IsDofManager >;

public:  
  //! type of Grid this DofManager belongs to 
  typedef Grid GridType;

public:
  typedef typename GridObjectStreamOrDefault<
    GridType, DummyObjectStream>::ObjectStreamType ObjectStreamType;

  typedef DataCollectorInterface<GridType, ObjectStreamType> DataCollectorType;

private:  
  typedef std::list< ManagedDofStorageInterface* > ListType;
  typedef typename ListType::iterator ListIteratorType;
  typedef typename ListType::const_iterator ConstListIteratorType;

  typedef LocalInterface< int > MemObjectCheckType;
  
  typedef std::list< ManagedIndexSetInterface * > IndexListType;
  typedef typename IndexListType::iterator IndexListIteratorType;
  typedef typename IndexListType::const_iterator ConstIndexListIteratorType;

  // list with MemObjects, for each DiscreteFunction we have one MemObject
  ListType memList_;

  // list of all different indexsets 
  IndexListType indexList_;

  // the dofmanager belong to one grid only 
  const GridType &grid_;

  // index set for mapping 
  mutable DataCollectorType dataInliner_;
  mutable DataCollectorType dataXtractor_;

  //! type of IndexSet change interfaces 
  typedef LocalInterface<typename GridType::
        template Codim<0>::Entity> LocalIndexSetObjectsType;

  mutable LocalIndexSetObjectsType indexSets_; 

  mutable LocalIndexSetObjectsType insertIndices_;
  mutable LocalIndexSetObjectsType removeIndices_;

  // lists containing all MemObjects 
  // to have fast access during resize and reserve 
  mutable MemObjectCheckType resizeMemObjs_;
  mutable MemObjectCheckType reserveMemObjs_;

  //! if chunk size if small then defaultChunkSize is used 
  const int defaultChunkSize_; 

  //! number of sequence, incremented every resize is called
  int sequence_; 
  
public: 
  typedef IndexSetRestrictProlong< ThisType , LocalIndexSetObjectsType >
    IndexSetRestrictProlongType;
  // this class needs to call resizeMemory 
  friend class IndexSetRestrictProlong< ThisType , LocalIndexSetObjectsType > ;

private:
  // combine object holding all index set for restrict and prolong 
  IndexSetRestrictProlongType indexRPop_; 
  
  //! memory over estimation factor for re-allocation 
  double memoryFactor_;
  //**********************************************************
  //**********************************************************
  //! Constructor 
  inline explicit DofManager ( const GridType *grid ) 
  : grid_( *grid ),
    defaultChunkSize_( 128 ),
    sequence_( 0 ),
    indexRPop_( *this, insertIndices_ , removeIndices_ ),
    memoryFactor_( Parameter :: getValidValue
      ( "fem.dofmanager.memoryfactor",  double( 1.1 ),
        ValidateNotLess< double >( 1.0 ) ) )
  {
    if( Parameter :: verbose() && (grid_.comm().rank() == 0) )
      std :: cout << "DofManager: Created for " << grid_.name()
                  << " with memory factor " << memoryFactor_
                  << "." << std :: endl;
  }

  // copy of dofmanagers is forbidden 
  DofManager( const ThisType & )
  {
    std::cerr << "DofManager(const DofManager &) not allowed!" << std :: endl;
    abort();
  }
  
  //! Desctructor, removes all MemObjects and IndexSetObjects 
  ~DofManager (); 

public:
  //! return factor to over estimate new memory allocation 
  double memoryFactor() const { return memoryFactor_; }
  
  //! add new index set to the list of the indexsets of this dofmanager
  template <class IndexSetType>
  inline void addIndexSet (const IndexSetType &iset); 

  //! remove index set from the indexsets list of this dofmanager
  template <class IndexSetType>
  inline void removeIndexSet (const IndexSetType &iset); 

  /** \brief add a managed dof storage to the dof manager. 
      \param dofStorage  dof storage to add which must fulfill the 
             ManagedDofStorageInterface 
  */
  template <class ManagedDofStorageImp>
  void addDofStorage(ManagedDofStorageImp& dofStorage);

  /** \brief remove a managed dof storage from the dof manager. 
      \param dofStorage  dof storage to remove which must fulfill the 
             ManagedDofStorageInterface 
  */
  template <class ManagedDofStorageImp>
  void removeDofStorage(ManagedDofStorageImp& dofStorage);

  //! returns the index set restrinction and prolongation operator
  IndexSetRestrictProlongType & indexSetRPop () 
  {
    // hier muss statt dessen ein Combiniertes Object erzeugt werden. 
    // dafuer sollte bei einhaengen der IndexSets ein Methoden Pointer
    // erzeugt werden, welcher die den IndexSet mit einem anderen Object
    // kombiniert 
    return indexRPop_;
  }

  //! if dofmanagers list is not empty return true 
  bool hasIndexSets() const 
  {
    return ! insertIndices_.empty(); 
  }
   
  /** \brief return used memory size of all MemObjects in bytes. */
  size_t usedMemorySize () const 
  {
    size_t used = 0;
    ConstListIteratorType endit = memList_.end();
    for(ConstListIteratorType it = memList_.begin(); it != endit ; ++it)
    {
      used += (*it)->usedMemorySize(); 
    }
    return used;
  }
  
  /** \brief resize memory before data restriction 
      during grid adaptation is done.
  */ 
  void resizeForRestrict () 
  {
    ++sequence_;
    resizeMemory();
  }
  
  /** \brief reserve memory for at least nsize elements 
      this will increase the sequence counter by 1 
      if useNsize is true, then nsize will be used as chunk size 
      otherwise max( nsize, defaultChunkSize_ )
  */ 
  void reserveMemory (int nsize, bool useNsize = false ) 
  {
    ++sequence_;
    int localChunkSize = (useNsize) ? nsize : std::max(nsize, defaultChunkSize_ );
    assert( localChunkSize > 0 );

    // reserves (size + chunkSize * elementMemory), see above 
    reserveMemObjs_.apply ( localChunkSize );
  }

  /** \brief return number of sequence, if dofmanagers memory was changed by
      calling some method like resize, then also this number will increase
     \note The increase of this number could be larger than 1 
  */
  int sequence () const { return sequence_; }

  /** \brief Resize index sets and memory due to what the mapper has as new size.
      \note This will increase the sequence counter by 1. 
  */
  void resize()
  {
    // new number in grid series 
    ++sequence_;

    IndexListIteratorType endit = indexList_.end();
    for(IndexListIteratorType it = indexList_.begin(); it != endit; ++it)
    {
      (*it)->resize(); 
    }
    resizeMemory();
  }

  /** \brief Inserts entity to all index sets added to dof manager. */
  template <class EntityType>
  inline void insertEntity(EntityType & en )
  {
    // insert new index 
    insertIndices_.apply( en );

    // resize memory 
    resizeMemory();
  }
          
  /** \brief Removes entity from all index sets added to dof manager. */
  template <class EntityType>
  inline void removeEntity(EntityType & en )
  {
    removeIndices_.apply( en );
  }

protected:  
  //! resize the MemObject if necessary 
  void resizeMemory()
  {
    int dummy = -1;
    // pass dummy parameter 
    resizeMemObjs_.apply ( dummy ); 
  }
  
public:
  /** \brief Compress all data that is hold by this dofmanager 
      \note This will increase the sequence counter by 1.
      \note This method is deprecated! 
  */
  void dofCompress() 
  {
    compress();
  }
  
  /** \brief Compress all data that is hold by this dofmanager 
      \note  This will increase the sequence counter by 1.
  */
  void compress() 
  {
    // mark next sequence 
    ++sequence_;

    // compress indexsets first 
    {
      IndexListIteratorType endit  = indexList_.end();
      for(IndexListIteratorType it = indexList_.begin(); it != endit; ++it)
      {
        // reset compressed so the next time compress of index set is called 
        (*it)->compress(); 
      }
    }

    // compress all data now 
    {
      ListIteratorType endit  = memList_.end();
      for(ListIteratorType it = memList_.begin(); it != endit ; ++it)
      {
        // if correponding index was not compressed yet, this is called in
        // the MemObject dofCompress, if index has not changes, nothing happens  
        // if IndexSet actual needs  no compress, nothing happens to the
        // data either 
        // also data is resized, which means the vector is getting shorter
        (*it)->dofCompress () ;
      }
    }
  }

  //! add data handler for data inlining to dof manager
  template <class DataCollType>
  void addDataInliner ( DataCollType & d)
  {
    dataInliner_ += d;
  }

  //! clear data inliner list 
  void clearDataInliners ()
  {
    dataInliner_.clear();
  }

  //! add data handler for data xtracting to dof manager
  template <class DataCollType>
  void addDataXtractor ( DataCollType & d)
  {
    dataXtractor_ += d;
  }

  //! clear data xtractor list 
  void clearDataXtractors ()
  {
    dataXtractor_.clear();
  }

  //! packs all data of this entity en and all child entities  
  template <class ObjectStreamType, class EntityType>
  void inlineData ( ObjectStreamType & str, EntityType & en )
  {
    dataInliner_.apply(str,en);
  }

  //! unpacks all data of this entity from message buffer 
  template <class ObjectStreamType, class EntityType>
  void xtractData ( ObjectStreamType & str, EntityType & en, size_t newElements )
  {
    // reserve memory for new elements 
    reserveMemory(newElements , true );
    // here the elements already have been created 
    // that means we can xtract data
    dataXtractor_.apply(str,en);
  }

private:
  //! only called from DofManagerFactory 
  //********************************************************
  // read-write Interface for index set 
  //********************************************************

  //! writes all underlying index sets to a file 
  bool writeIndexSets(const std::string filename, int timestep);
  //! reads all underlying index sets from a file 
  bool readIndexSets(const std::string filename, int timestep);

  // generate index set filename 
  std::string generateIndexSetName(const std::string& filename,
                                   const int count) const
  {
    std::string newFilename (filename);
    newFilename += "_";
    // add number 
    {
      std::stringstream tmp;
      tmp << count;
      newFilename += tmp.str();
    }

    newFilename += "_";
    return newFilename;
  }
}; // end class DofManager

//***************************************************************************
//
//  inline implemenations 
//
//***************************************************************************

template <class GridType>
inline DofManager<GridType>::~DofManager () 
{
  if(memList_.size() > 0)
  {
    while( memList_.rbegin() != memList_.rend())
    {
      DofStorageInterface * mobj = (* memList_.rbegin() );
      memList_.pop_back();
      
      // alloc new mem an copy old mem 
      dverb << "Removing '" << mobj->name() << "' from DofManager!\n";  
    }
  }

  if(indexList_.size() > 0)
  {
    while ( indexList_.rbegin() != indexList_.rend()) 
    {
      ManagedIndexSetInterface* iobj = (* indexList_.rbegin() );
      indexList_.pop_back();
      if(iobj) delete iobj;
    }
  }
}

template <class GridType>
template <class IndexSetType>
inline void DofManager<GridType>::
addIndexSet (const IndexSetType &iset)
{
  typedef typename GridType::template Codim<0>::Entity EntityType;
  typedef ManagedIndexSet< IndexSetType, EntityType > ManagedIndexSetType;
  
  ManagedIndexSetType * indexSet = 0;
  
  IndexListIteratorType endit = indexList_.end();
  for(IndexListIteratorType it = indexList_.begin(); it != endit; ++it)
  {
    // check equality 
    // and increase counter if equal
    if( (*it)->increaseReference(iset) )
    {
      indexSet = static_cast<ManagedIndexSetType *> ((*it));
      break;
    }
  }
  
  if(!indexSet) 
  { 
    indexSet = new ManagedIndexSetType ( iset, indexSets_ , insertIndices_ , removeIndices_  );
      
    ManagedIndexSetInterface* iobj = indexSet;
    // push to front to search latest index sets fast
    indexList_.push_front( iobj );
  }
  return ; 
}

template <class GridType>
template <class IndexSetType>
inline void DofManager<GridType>::
removeIndexSet (const IndexSetType &set)
{
  // search object in list an remove it
  IndexListIteratorType endit = indexList_.end();
  for(IndexListIteratorType it = indexList_.begin(); it != endit; ++it)
  {
    // decrease reference counter 
    // and delete if refernce is zero
    if( (*it)->decreaseReference( set ) )
    {
      // get obj pointer  
      ManagedIndexSetInterface* set = *it;
      // remove from list 
      indexList_.erase( it );
      // delete proxy 
      delete set;
      return ;
    }
  }

  // we should never get here
  DUNE_THROW(InvalidStateException,"Could not remove index set!");
}

template <class GridType>
template <class ManagedDofStorageImp>
void
DofManager<GridType>::
addDofStorage(ManagedDofStorageImp& dofStorage)
{
  dverb << "Adding '" << dofStorage.name() << "' to DofManager! \n";

  // make sure we got an ManagedDofStorage
  ManagedDofStorageInterface* obj = &dofStorage;

  // push_front, makes search faster 
  memList_.push_front( obj );

  // add the special object to the memResize list object 
  resizeMemObjs_  += dofStorage.resizeMemoryObject();

  // the same for the reserve call  
  reserveMemObjs_ += dofStorage.reserveMemoryObject();
}


template <class GridType>
template <class ManagedDofStorageImp>
void
DofManager<GridType>::
removeDofStorage(ManagedDofStorageImp& dofStorage)
{
  // make sure we got an ManagedDofStorage
  ManagedDofStorageInterface* obj = &dofStorage;

  // search list starting from tail  
  ListIteratorType endit = memList_.end();
  for( ListIteratorType it = memList_.begin();
       it != endit ; ++it)
  {
    if(*it == obj)
    {
      // alloc new mem and copy old mem 
      memList_.erase( it );

      dvverb << "Remove '" << dofStorage.name() << "' from DofManager!\n";

      // remove from list 
      resizeMemObjs_.remove( dofStorage.resizeMemoryObject() );
      reserveMemObjs_.remove( dofStorage.reserveMemoryObject() );

      return ;
    }
  }
}

template <class GridType>
inline bool DofManager<GridType>::
writeIndexSets(const std::string filename , int timestep )
{
  int count = 0;
  IndexListIteratorType endit = indexList_.end();
  for(IndexListIteratorType it = indexList_.begin(); it != endit; ++it)
  {
    std::string newFilename = generateIndexSetName(filename, count);
    (*it)->write_xdr(newFilename.c_str(),timestep); 
    ++ count;
  }
  return true;
}

template <class GridType>
inline bool DofManager<GridType>::
readIndexSets(const std::string filename , int timestep )
{
  int count = 0;
  IndexListIteratorType endit = indexList_.end();
  for(IndexListIteratorType it = indexList_.begin(); it != endit; ++it)
  {
    // create filename 
    std::string newFilename = generateIndexSetName(filename, count);
    
    // create filename that is used by index sets 
    std::string fileToRead = genFilename("", newFilename.c_str(), timestep);
      
    // check if file exists, and skip if not 
    FILE * testfile = fopen(fileToRead.c_str(),"r");
    if( testfile )
    {
      fclose( testfile );
      (*it)->read_xdr(newFilename.c_str(),timestep);
    }
    else if( Parameter :: verbose() )
    {
      std::cout << "WARNING: Skipping " << newFilename << " in DofManager::read_xdr!" << std::endl;
    }
    ++count;
  }
  return true;
}

//@} 



  /** \class DofManagerFactory
   *  \ingroup DofManager
   *  \brief Singleton provider for the DofManager
   *
   *  DofManagerFactory guarantees that at most one instance of DofManager
   *  is generated for each grid.
   */
  template< class DofManager >
  class DofManagerFactory
  {
    typedef DofManagerFactory< DofManager > ThisType;

  public:
    typedef DofManager DofManagerType;
    typedef typename DofManagerType :: GridType GridType; 

  private:
    typedef const GridType *KeyType;

    typedef SingletonList< KeyType, DofManagerType > DMProviderType;

  public:
    /** \brief obtain a reference to the DofManager for a given grid
     *
     *  \param[in]  grid  grid for which the DofManager is desired
     *
     *  \returns a reference to the singleton instance of the DofManager
     */
    inline static DofManagerType &getDofManager ( const GridType &grid )
    {
      DofManagerType *dm = getDmFromList( grid );
      if( !dm )
        return DMProviderType :: getObject( &grid );
      return *dm;
    } 

    //! delete the dof manager that belong to the given grid 
    inline static void deleteDofManager ( DofManagerType &dm )
    {
      DMProviderType :: removeObject( &dm );
    }

    //! writes DofManager of corresponding grid, when DofManager exists 
    inline static bool 
    writeDofManager ( const GridType &grid,
                      const std :: string &filename,
                      int timestep )
    {
      DofManagerType *dm = getDmFromList( grid );
      if( dm )
        return dm->writeIndexSets( filename, timestep );
      return false;
    }

    //! reads DofManager of corresponding grid, when DofManager exists 
    inline static bool 
    readDofManager ( const GridType &grid,
                     const std :: string &filename,
                     int timestep )
    {
      DofManagerType *dm = getDmFromList( grid );
      if( dm )
        return dm->readIndexSets( filename, timestep );
      return false;
    }

  private: 
    // return pointer to dof manager for given grid 
    inline static DofManagerType *getDmFromList( const GridType &grid )
    {
      return (DMProviderType :: getObjFromList( &grid )).first;
    }
  };


} // end namespace Dune 

#endif
