#ifndef DUNE_DOFMAPPERINTERFACE_HH
#define DUNE_DOFMAPPERINTERFACE_HH

//- Dune includes
#include <dune/fem/misc/bartonnackmaninterface.hh>
#include <dune/fem/space/common/basefunctioninterface.hh>

namespace Dune {

/** @addtogroup DofMapper  

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
template< class DofMapperTraits >
class DofMapperInterface
: public BartonNackmanInterface< DofMapperInterface< DofMapperTraits >,
                                 typename DofMapperTraits :: DofMapperType >
{
public:
  typedef DofMapperTraits Traits;

  typedef typename Traits :: DofMapperType DofMapperType;
  
private:
  typedef DofMapperInterface< Traits > ThisType;
  typedef BartonNackmanInterface< ThisType, DofMapperType > BaseType;

public:
  typedef typename Traits :: EntityType EntityType;

  typedef typename Traits :: DofMapIteratorType DofMapIteratorType;

protected:
  using BaseType :: asImp;

public: 
  //! return number of dofs for special function space and grid on
  //! specified level
  int size () const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().size());
    return asImp().size();
  }

  /** \brief obtain a begin iterator for the mapping on an entity
   *
   *  To obtain the entire mapping of local DoF numbers to global DoF numbers,
   *  using these iterators can be much more efficient than calling
   *  mapToGlobal for each local DoF.
   *
   *  \param[in]  entity  entity (codimension 0) for which the mapping shall
   *                      be iterated
   *
   *  \returns begin iterator for the local DoF mapping
   */
  inline DofMapIteratorType begin ( const EntityType &entity ) const
  {
    CHECK_INTERFACE_IMPLEMENTATION( asImp().begin( entity ) );
    return asImp().begin( entity );
  }
  
  /** \brief obtain an end iterator for the mapping on an entity
   *
   *  To obtain the entire mapping of local DoF numbers to global DoF numbers,
   *  using these iterators can be much more efficient than calling
   *  mapToGlobal for each local DoF.
   *
   *  \param[in]  entity  entity (codimension 0) for which the mapping shall
   *                      be iterated
   *
   *  \returns end iterator for the local DoF mapping
   */
  inline DofMapIteratorType end ( const EntityType &entity ) const
  {
    CHECK_INTERFACE_IMPLEMENTATION( asImp().end( entity ) );
    return asImp().end( entity );
  }

  /** return new size of space, i.e. after adaptation 
    \returns new size of space 
  */
  int newSize() const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().newSize());
    return asImp().newSize();
  }

  /** \brief map a local DoF number to a global one
   *
   *  \param[in]  entity    entity the DoF belongs to
   *  \param[in]  localDof  local number of the DoF
   *
   *  \returns global number of the DoF
   */
  int mapToGlobal ( const EntityType &entity, const int localDof ) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().mapToGlobal(entity ,localDof));
    return asImp().mapToGlobal( entity , localDof );
  }
 
  /** \brief obtain maximal number of DoFs on one entity
   * 
   *  \note The naming of this method is misleading since the number of DoFs on
   *        one entity need not be constant (e.g. Lagrange base functions on a
   *        hybrid grid).
   */
  int numDofs () const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numDofs());
    return asImp().numDofs();
  }
  
  /** \brief map a local DoF number of an entity to a global one
   *
   *  \param[in]  entity    entity the DoF belongs to
   *  \param[in]  localDof  local number of the DoF
   *
   *  \returns global number of the DoF
   */
  template <class EntityImp> 
  int mapEntityDofsToGlobal ( const EntityImp &entity, const int localDof ) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(
        asImp().mapEntityDofsToGlobal(entity ,localDof));
    return asImp().mapEntityDofsToGlobal( entity , localDof );
  }
 
  /** \brief obtain number of DoFs on an entity
   * 
   *  \param[in]  entity  entity of codimension 0
   *  
   *  \returns number of DoFs on the entity
   */
  inline int numDofs ( const EntityType &entity ) const
  {
    CHECK_INTERFACE_IMPLEMENTATION( asImp().numDofs( entity ) );
    return asImp().numDofs( entity );
  }

  /** \brief return number of holes for data block */
  int numberOfHoles(const int block) const 
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numberOfHoles(block));
    return asImp().numberOfHoles(block); 
  }
  
  /** \brief return old index of hole for data block */
  int oldIndex (const int hole, const int block) const 
  { 
    CHECK_INTERFACE_IMPLEMENTATION(asImp().oldIndex(hole,block));
    return asImp().oldIndex(hole,block); 
  }
    
  /** \brief return new index of hole for data block */
  int newIndex (const int hole, const int block) const 
  { 
    CHECK_INTERFACE_IMPLEMENTATION(asImp().newIndex(hole,block));
    return asImp().newIndex(hole,block); 
  }

  /** \brief return true if compress will affect data */
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

  /** \brief return old offsets for given block */
  int oldOffSet(const int block) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().oldOffSet(block));
    return asImp().oldOffSet(block);
  }

  /** \brief return current offsets for given block */
  int offSet(const int block) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().offSet(block));
    return asImp().offSet(block);
  }

  /** return number of supported blocks */
  int numBlocks() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numBlocks());
    return asImp().numBlocks();
  }
};



template< class EntityImp, class DofMapperImp >
class DefaultDofMapIterator
{
public:
  typedef EntityImp EntityType;
  typedef DofMapperImp DofMapperType;

  enum IteratorType { beginIterator, endIterator };

private:
  typedef DefaultDofMapIterator< EntityType, DofMapperType > ThisType;

protected:
  const EntityType &entity_;
  const DofMapperType &dofMapper_;
  int dof_;
  
public:
  inline DefaultDofMapIterator ( const IteratorType type,
                                 const EntityType &entity,
                                 const DofMapperType &dofMapper )
  : entity_( entity ),
    dofMapper_( dofMapper ),
    dof_( type == beginIterator ? 0 : dofMapper_.numDofs( entity ) )
  {}

  inline DefaultDofMapIterator ( const ThisType &other )
  : entity_( other.entity_ ),
    dofMapper_( other.dofMapper_ ),
    dof_( other.dof_ )
  {}

  inline ThisType &operator++ ()
  {
    ++dof_;
    return *this;
  }

  inline bool operator== ( const ThisType &other ) const
  {
    return dof_ == other.dof_;
  }

  inline bool operator!= ( const ThisType &other ) const
  {
    return dof_ != other.dof_;
  }

  inline int local () const
  {
    return dof_;
  }

  inline int global () const
  {
    return dofMapper_.mapToGlobal( entity_, dof_ );
  }
};



//! Default implementation for DofMappers, empty at this moment
template< class DofMapperTraits >
class DofMapperDefault
: public DofMapperInterface< DofMapperTraits >
{
public:
  typedef DofMapperTraits Traits;

  typedef typename Traits :: EntityType EntityType;

private:
  typedef DofMapperDefault< Traits > ThisType;
  typedef DofMapperInterface< Traits > BaseType;

protected:
  using BaseType :: asImp;
  
public:
  /** \copydoc DofMapperInterface::numDofs( const EntityType &entity ) const
      \note This implementation just returns number of all dofs 
  */
  inline int numDofs ( const EntityType &entity ) const
  {
    return asImp().numDofs ();
  }

  //! update mapper, default does nothing 
  void update () {}

  //! return old offsets for block number, default returns zero 
  int oldOffSet(const int block) const { return 0; }

  //! return current offsets for block number, default returns zero 
  int offSet(const int block) const { return 0; }

  //! return number of supported blocks, default is 1  
  int numBlocks() const { return 1; }
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
