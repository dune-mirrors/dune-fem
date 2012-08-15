#ifndef DUNE_FEM_DOFMAPPER_HH
#define DUNE_FEM_DOFMAPPER_HH

//- Dune includes
#include <dune/fem/misc/bartonnackmaninterface.hh>
#include <dune/fem/gridpart/emptyindexset.hh>

namespace Dune
{

  /** @addtogroup DofMapper  

      Dof Mapper are special implementations to provide the mapToGlobal
      feature of the spaces. Furthermore, during adaptation the mapper
      provide information about holes in the data vectors. 

      \remarks 
      The DofMapper interface is described by the class
      DofMapper.

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
  class DofMapper
  : public Fem::BartonNackmanInterface< DofMapper< DofMapperTraits >,
                                   typename DofMapperTraits :: DofMapperType >
  {
    typedef DofMapper< DofMapperTraits > ThisType;
    typedef Fem::BartonNackmanInterface< ThisType, typename DofMapperTraits::DofMapperType > BaseType;

  public:
    typedef DofMapperTraits Traits;

    //! type of the DofMapper implementation
    typedef typename Traits::DofMapperType DofMapperType;
   
    //! type of codimension 0 entities
    typedef typename Traits::ElementType ElementType;

    typedef ElementType EntityType;

    //! type of the dof map iterator
    typedef typename Traits::DofMapIteratorType DofMapIteratorType;

  protected:
    using BaseType::asImp;

  public: 
    /** \brief  return number of dofs for special function space and grid on
                specified level */
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
    DofMapIteratorType begin ( const ElementType &entity ) const
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
    DofMapIteratorType end ( const ElementType &entity ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().end( entity ) );
      return asImp().end( entity );
    }

    /** \brief returns true if DoFs for given codimension exist 
     *
     *  \param[in]  codim   codimension to check 
     *
     *  \returns true if DoFs for codimension exist 
     */
    bool contains ( const int codim ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().contains( codim ) );
      return asImp().contains( codim );
    }

    /** \brief map each local DoF number to a global one
     *
     *  \param[in]  element  element, the DoFs belong to
     *  \param[in]  f        functor to call for each DoF
     *
     *  The functor has to be a copyable object satisfying the following
     *  interface:
     *  \code
     *  struct Functor
     *  {
     *    // application operator
     *    void operator() ( int localDoF, int globalDoF );
     *  };
     *  \endcode
     *
     *  For each DoF to be mapped, this method will call the application operator
     *  once.
     *  
     *  \note There is no guarantee on the order, in which the functor is applied.
     */
    template< class Functor >
    void mapEach ( const ElementType &element, Functor f ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().mapEach( element, f ) );
    }

    /** \brief map a local DoF number to a global one
     *
     *  \param[in]  entity    entity the DoF belongs to
     *  \param[in]  localDof  local number of the DoF
     *
     *  \returns global number of the DoF
     */
    int mapToGlobal ( const ElementType &entity, const int localDof ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().mapToGlobal(entity ,localDof));
      return asImp().mapToGlobal( entity , localDof );
    }
    
     /** \brief map a local DoF number of an entity to a global one
     *
     *  \param[in]  entity    entity the DoF belongs to
     *  \param[in]  localDof  local number of the DoF
     *
     *  \returns global number of the DoF
     */
    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION
        ( asImp().mapEntityDofToGlobal(entity ,localDof) );
      return asImp().mapEntityDofToGlobal( entity , localDof );
    }
    
    /** \brief obtain maximal number of DoFs on one entity
     */
    int maxNumDofs () const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().maxNumDofs());
      return asImp().maxNumDofs();
    }

    /** \brief obtain number of DoFs on an entity
     * 
     *  \param[in]  element  entity of codimension 0
     *  
     *  \returns number of DoFs on the entity
     */
    int numDofs ( const ElementType &element ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().numDofs( element ) );
      return asImp().numDofs( element );
    }

    /** \brief obtain number of DoFs actually belonging to an entity
     *
     *  In contrast to numDofs, this method returns the number of DoFs actually
     *  associated with an entity (usually a subentity). We have the following
     *  relation for an entity \f$E\f$ of codimension 0:
     *  \f[
     *  \mathrm{numDofs}( E ) = \sum_{e \subset E} \mathrm{numEntityDofs}( e ),
     *  \f]
     *  where \f$\subset\f$ denotes the subentity relation.
     * 
     *  \param[in]  entity  entity of codimension
     *  
     *  \returns number of DoFs on the entity
     */
    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().numEntityDofs( entity ) );
      return asImp().numEntityDofs( entity );
    }

    /** \brief Check, whether the data in a codimension has fixed size */
    bool fixedDataSize ( int codim ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().fixedDataSize( codim ) );
      return asImp().fixedDataSize( codim );
    }

    /** \brief return number of holes for data block */
    int numberOfHoles(const int block) const 
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().numberOfHoles(block));
      return asImp().numberOfHoles(block); 
    }
    
    /** \brief return old index of hole for data block (with resprect to new offset) */
    int oldIndex (const int hole, const int block) const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().oldIndex(hole,block));
      return asImp().oldIndex(hole,block); 
    }
      
    /** \brief return new index of hole for data block (with resprect to new offset) */
    int newIndex (const int hole, const int block) const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().newIndex(hole,block));
      return asImp().newIndex(hole,block); 
    }

    /** \brief return true if compress will affect data */
    bool consecutive () const 
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().consecutive() );
      return asImp().consecutive();
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

    /** \brief return number of supported blocks */
    int numBlocks() const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().numBlocks());
      return asImp().numBlocks();
    }
  };



  template< class Element, class DofMapper >
  class DefaultDofMapIterator
  {
    typedef DefaultDofMapIterator< Element, DofMapper > ThisType;

  public:
    typedef Element ElementType;
    typedef DofMapper DofMapperType;

    enum IteratorType { beginIterator, endIterator };
    
    DefaultDofMapIterator ( const IteratorType type,
                            const ElementType &entity,
                            const DofMapperType &dofMapper )
    : entity_( entity ),
      dofMapper_( dofMapper ),
      dof_( type == beginIterator ? 0 : dofMapper_.numDofs( entity ) )
    {}

    DefaultDofMapIterator ( const ThisType &other )
    : entity_( other.entity_ ),
      dofMapper_( other.dofMapper_ ),
      dof_( other.dof_ )
    {}

    const ThisType &operator++ ()
    {
      ++dof_;
      return *this;
    }

    bool operator== ( const ThisType &other ) const
    {
      return dof_ == other.dof_;
    }

    bool operator!= ( const ThisType &other ) const
    {
      return dof_ != other.dof_;
    }

    int local () const
    {
      return dof_;
    }

    int global () const
    {
      return dofMapper_.mapToGlobal( entity_, dof_ );
    }

  protected:
    const ElementType &entity_;
    const DofMapperType &dofMapper_;
    int dof_;
  };



  //! Default implementation for DofMappers, empty at this moment
  template< class DofMapperTraits >
  class DofMapperDefault
  : public DofMapper< DofMapperTraits >
  {
    typedef DofMapperDefault< DofMapperTraits > ThisType;
    typedef DofMapper< DofMapperTraits > BaseType;

  public:
    typedef DofMapperTraits Traits;

    typedef typename Traits::ElementType ElementType;

  protected:
    using BaseType :: asImp;

    //! for dune-fem index set use new method consecutive
    template< class IndexSetType >
    bool checkConsecutive ( const IndexSetType &indexSet ) const
    {
      return indexSet.consecutive();
    }

  public:
    /** \brief Check, whether the data in a codimension has fixed size */
    bool fixedDataSize ( int codim ) const { return false; }

    //! return old offsets for block number, default returns zero 
    int oldOffSet(const int block) const { return 0; }

    //! return current offsets for block number, default returns zero 
    int offSet(const int block) const { return 0; }

    //! return number of supported blocks, default is 1  
    int numBlocks() const { return 1; }
  }; 

  /// @}

} // end namespace Dune

#endif // #ifndef DUNE_FEM_DOFMAPPER_HH
