#ifndef DUNEFEM_INDEXSETS_HH
#define DUNEFEM_INDEXSETS_HH

//- system includes 
#include <iostream>
#include <string> 
#include <rpc/xdr.h>
#include <cassert>

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/grid/common/indexidset.hh>

//- Dune fem includes 
#include <dune/fem/gridpart/emptyindexset.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/io/file/persistencemanager.hh>

/** @file
 @brief Provides default index set class for persistent index sets. 
*/

namespace Dune
{

  /** \brief  class implementing the DUNE grid index set interface for any DUNE
              fem index set. However, one should not cast to this class as
              an interface class.
  */
  template< class GridImp, class Imp >
  class DuneGridIndexSetAdapter 
  : public BartonNackmanInterface< DuneGridIndexSetAdapter< GridImp, Imp >, Imp >,
    public EmptyIndexSet,
    public IndexSet< GridImp, Imp >
  {
    typedef DuneGridIndexSetAdapter< GridImp, Imp > ThisType;
    typedef BartonNackmanInterface< ThisType, Imp > BaseType;

    friend class Conversion< ThisType, EmptyIndexSet >;

  protected:
    using BaseType::asImp;

    typedef IndexSet< GridImp, Imp > DuneIndexSetType;

  public:
    //! type of grid 
    typedef GridImp GridType;

    //! type of index (i.e. unsigned int)
    typedef typename DuneIndexSetType::IndexType IndexType;

    using DuneIndexSetType::index; 

    //! type of codimension 0 entity 
    typedef typename GridType::template Codim< 0 >::Entity EntityCodim0Type; 

    //! constructor storing grid reference 
    explicit DuneGridIndexSetAdapter ( const GridType &grid )
    : grid_(grid) 
    {}

    //! copy constructor 
    DuneGridIndexSetAdapter ( const ThisType &other )
    : grid_( other.grid_ )
    {}
    
  public:
    //****************************************************************
    //
    //  INTERFACE METHODS for DUNE INDEX SETS 
    //
    //****************************************************************
    //! return global index of entity 
    template< class EntityType >
    IndexType index ( const EntityType &entity ) const
    {
      // return index of entity 
      enum { codim = EntityType::codimension };
      return this->template index< codim >( entity, 0 );
    }

    //! return subIndex of given entity
    template< int codim >
    IndexType subIndex ( const EntityCodim0Type &entity, const int localNum ) const
    {
      // return sub index of entity 
      return this->template index< codim >( entity, localNum );
    }

    //////////////////////////////////////////////////////////////////
    //
    //  DUNE fem index method implementer interface  
    //
    //////////////////////////////////////////////////////////////////
  protected:
    //! return index for entity  
    template< int codim, class EntityType >
    IndexType indexImp ( const EntityType &entity, const int localNum ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().template indexImp< codim >( entity, localNum ) );
      return asImp().template indexImp< codim >( entity, localNum );
    } 

  public:  
    //////////////////////////////////////////////////////////////////
    //
    //  DUNE fem index method user interface  
    //
    //////////////////////////////////////////////////////////////////
    //! return index for entity  
    template< int codim, class EntityType >
    IndexType index ( const EntityType &entity, const int localNum ) const
    {
      return this->template indexImp< codim >( entity, localNum );
    } 

    //! insert new index for entity to set 
    void insertEntity ( const EntityCodim0Type &entity )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().insertEntity( entity ) );
    }

    //! remove index for entity from index set 
    void removeEntity ( const EntityCodim0Type &entity )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().removeEntity( entity ) );
    }
    
  protected:  
    //! reference to grid 
    const GridType &grid_;
  };

  /** \brief ConsecutivePersistentIndexSet is the base class for all index sets that 
      are consecutive and also persistent. Implementations of this type
      are for example AdaptiveLeafIndexSet and DGAdaptiveLeafIndexSet. 
  */
  template< class GridImp, class Imp >
  class ConsecutiveIndexSet
  : public DuneGridIndexSetAdapter< GridImp, Imp >
  {
    typedef ConsecutiveIndexSet< GridImp, Imp > ThisType;
    typedef DuneGridIndexSetAdapter< GridImp, Imp > BaseType;

    typedef Imp ImplementationType;

  public:  
    //! type of grid 
    typedef GridImp GridType;

  protected:
    //! Conschdrugdor 
    explicit ConsecutiveIndexSet ( const GridType &grid )
    : BaseType( grid )
    {}

  private:
    // no copying & no assignment
    ConsecutiveIndexSet ( const ThisType & );
    ThisType &operator= ( const ThisType & );

  public:
    //! returns true since we deal with a consecutive index set 
    bool consecutive () const
    {
      return true;
    }

    //! return true if the index set is persistent 
    bool persistent () const
    {
      return false;
    }

    //! remove holes and make index set consecutive 
    bool compress()
    {
      return asImp().compress();
    } 

  protected:
    // use asImp from BaseType
    using BaseType::asImp;
  };



  /** \brief ConsecutivePersistentIndexSet is the base class for 
      all index sets that are persistent. Implementations of this type
      are for example all ConsecutivePersistenIndexSets. 
  */
  template< class GridImp, class Imp >
  class PersistentIndexSet
  : public DuneGridIndexSetAdapter< GridImp, Imp > ,
    public PersistentObject 
  {
    typedef PersistentIndexSet< GridImp, Imp > ThisType;
    typedef DuneGridIndexSetAdapter< GridImp, Imp > BaseType;

    typedef Imp ImplementationType;

  public:  
    //! type of entity with codimension 0
    typedef typename BaseType::EntityCodim0Type EntityCodim0Type;
    
    //! type of grid 
    typedef GridImp GridType;
    //! type of DoF manager
    typedef DofManager< GridType > DofManagerType;

  protected:
    /** \brief constructor */
    explicit PersistentIndexSet ( const GridType &grid )
      // here false, because methods have to be overloaded
      : BaseType(grid),
        dofManager_( DofManagerType::instance( grid ) )
    {
      // add persistent index set to dofmanagers list 
      dofManager_.addIndexSet( asImp() );
    }

  private:
    // no copying & no assignment
    PersistentIndexSet ( const ThisType & );
    ThisType &operator= ( const ThisType & );

  public:
    //! destructor remoing index set from dof manager  
    ~PersistentIndexSet () 
    {
      // remove persistent index set from dofmanagers list 
      dofManager_.removeIndexSet( asImp() );
    }

    //! return true if the index set is persistent 
    bool persistent () const
    {
      return true;
    }

  protected:
    friend class PersistenceManager ;

    /** \copydoc Dune::PersistentObject :: backup */
    virtual void backup() const 
    {
      std::string filename ( PersistenceManager :: uniqueFileName ( "indexset" ) );
#ifndef NDEBUG 
      bool success = 
#endif
        write_xdr( filename );
      assert( success );
    }

    /** \copydoc Dune::PersistentObject :: restore */
    virtual void restore() 
    {
      std::string filename ( PersistenceManager :: uniqueFileName ( "indexset" ) );
#ifndef NDEBUG 
      bool success = 
#endif
        read_xdr( filename );
      assert( success );
    }

    //! write index set to file 
    virtual bool write_xdr( const std::string & ) const { return false; }
    //! read index set from file 
    virtual bool read_xdr( const std::string & ) { return false; }

    using BaseType::asImp;

    // reference to dof manager 
    DofManagerType& dofManager_;
  };


  /** \brief ConsecutivePersistentIndexSet is the base class for all index sets that 
      are consecutive and also persistent. Implementations of this type
      are for example AdaptiveLeafIndexSet and DGAdaptiveLeafIndexSet. 
  */
  template< class GridImp, class Imp >
  class ConsecutivePersistentIndexSet
  : public PersistentIndexSet< GridImp, Imp >
  {
    typedef ConsecutivePersistentIndexSet< GridImp, Imp > ThisType;
    typedef PersistentIndexSet< GridImp, Imp > BaseType;

    typedef Imp ImplementationType;

  public:  
    //! type of grid 
    typedef GridImp GridType;

  protected:
    //! Conschdrugdor 
    explicit ConsecutivePersistentIndexSet ( const GridType &grid )
    : BaseType( grid )
    {}

  private:
    // no copying & no assignment
    ConsecutivePersistentIndexSet ( const ThisType & );
    ThisType &operator= ( const ThisType & );

  public:
    //! returns true since we deal with a consecutive index set 
    bool consecutive () const
    {
      return true;
    }

    //! remove holes and make index set consecutive 
    bool compress()
    {
      return asImp().compress();
    } 

  protected:
    // use asImp from BaseType
    using BaseType::asImp;
  };

} // end namespace Dune

#endif
