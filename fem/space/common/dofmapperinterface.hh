#ifndef DUNE_DOFMAPPERINTERFACE_HH
#define DUNE_DOFMAPPERINTERFACE_HH

//- Dune includes
#include <dune/fem/space/common/basefunctioninterface.hh>

namespace Dune {

/** @defgroup DofMapper DofMapperInterface 
    @ingroup DiscreteFunctionSpace

    Dof Mapper are special implementations to provide the mapToGlobal
    feature of the spaces. Furthermore, during adaptation the mapper
    provide information about holes in the data vectors. 

    \remarks 
    The DofMapper interface is described by the class
    DofMapperInterface.

  @{
 */

//-------------------------------------------------------------------------
//-
//-  --MapperInterface 
//-
//-------------------------------------------------------------------------
/** \brief 
   Interface for calculating the size of a function space for a grid on a
   specified level.
   Furthermore the local to global mapping of dof number is done. 
   Also during grid adaptation this mapper knows about old and new indices
   of entities. 
*/
template <class DofMapperImp> 
class DofMapperInterface
{
public: 
  //! return number of dofs for special function space and grid on
  //! specified level
  int size () const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().size());
    return asImp().size();
  }

  //! map a local dof num of a given entity to a global dof num
  template <class EntityType>
  int mapToGlobal ( const EntityType &entity, const int localDof ) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().mapToGlobal(entity ,localDof));
    return asImp().mapToGlobal( entity , localDof );
  }

  //! return new size of space, i.e. after adaptation 
  int newSize() const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().newSize());
    return asImp().newSize();
  }
  
  //! return number of dofs on element
  int numDofs () const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numDofs());
    return asImp().numDofs();
  } 

  //! return number of holes for data block 
  int numberOfHoles(const int block) const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numberOfHoles(block));
    return asImp().numberOfHoles(block); 
  }
  
  //! return old index of hole for data block 
  int oldIndex (const int hole, const int block) const 
  { 
    CHECK_INTERFACE_IMPLEMENTATION(asImp().oldIndex(hole,block));
    return asImp().oldIndex(hole,block); 
  }
    
  //! return new index of hole for data block 
  int newIndex (const int hole, const int block) const 
  { 
    CHECK_INTERFACE_IMPLEMENTATION(asImp().newIndex(hole,block));
    return asImp().newIndex(hole,block); 
  }

  //! return true if compress will affect data  
  bool needsCompress () const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().needsCompress());
    return asImp().needsCompress ();
  }

  //! update mapper, 
  //! i.e. calculate new insertion points of blocks 
  void update ()
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().update());
  }

  //! return old offsets for given block 
  int oldOffSet(const int block) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().oldOffSet(block));
    return asImp().oldOffSet(block);
  }

  //! return current offsets for given block 
  int offSet(const int block) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().offSet(block));
    return asImp().offSet(block);
  }

  //! return number of supported blocks  
  int numBlocks() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numBlocks());
    return asImp().numBlocks();
  }
  
private:  
  //! Barton-Nackman trick 
  DofMapperImp &asImp()  { return static_cast<DofMapperImp &>(*this); };
  //! Barton-Nackman trick 
  const DofMapperImp &asImp() const { return static_cast<const DofMapperImp &>(*this); };
};

//! Default implementation for DofMappers, empty at this moment
template <class DofMapperImp> 
class DofMapperDefault : public DofMapperInterface<DofMapperImp>
{
public:
  //! update mapper, default does nothing 
  void update () {}

  //! return old offsets for block number, default returns zero 
  int oldOffSet(const int block) const { return 0; }

  //! return current offsets for block number, default returns zero 
  int offSet(const int block) const { return 0; }

  //! return number of supported blocks, default is 1  
  int numBlocks() const { return 1; }
  
private:  
  //! Barton-Nackman trick 
  DofMapperImp &asImp()  { return static_cast<DofMapperImp &>(*this); };
  //! Barton-Nackman trick 
  const DofMapperImp &asImp() const { return static_cast<const DofMapperImp &>(*this); };
}; 

//! Key for Mapper singleton list 
template <class IndexSetImp>
class MapperSingletonKey 
{
  const IndexSetImp & indexSet_; 
  const int numDofs_; 
public:
  //! constructor taking index set and numDofs 
  MapperSingletonKey(const IndexSetImp & indexSet, int numDofs ) : indexSet_(indexSet) ,  numDofs_(numDofs) {}
  //! copy constructor 
  MapperSingletonKey(const MapperSingletonKey &org) : indexSet_(org.indexSet_) , numDofs_(org.numDofs_) {}
  //! returns true if indexSet pointer and numDofs are equal 
  bool operator == (const MapperSingletonKey & otherKey) const 
  {
    return ((&indexSet_ == &otherKey.indexSet_) && (numDofs_ == otherKey.numDofs_));
  }

  //! return reference to index set 
  const IndexSetImp & indexSet() const { return indexSet_; }
  //! return number of dofs 
  const int numDofs () const { return numDofs_; }
};

//! Factory class for SingletonList to tell how objects are created and
//! how compared.
template <class KeyImp, class ObjectImp>
class MapperSingletonFactory
{
  public:
  //! create new mapper  
  static ObjectImp * createObject( const KeyImp & key )
  {
    // create Object of MapperType = ObjectImp 
    return new ObjectImp(key.indexSet(),key.numDofs());
  }
  //! delete mapper object 
  static void deleteObject( ObjectImp * obj )
  {
    delete obj;
  }
};

/// @}

} // end namespace Dune
#endif
