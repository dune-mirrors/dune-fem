#ifndef DUNE_DISCRETEFUNCTIONSPACE_HH
#define DUNE_DISCRETEFUNCTIONSPACE_HH

//- system includes
#include <assert.h>

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/fem/space/common/geometryconversion.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/basefunctioninterface.hh>
#include <dune/fem/function/common/temporarylocalfunction.hh>
#include <dune/fem/function/common/localfunctionwrapper.hh>


//- local includes 
#include "allgeomtypes.hh"
#include "singletonlist.hh"

namespace Dune
{

  /** @addtogroup DiscreteFunctionSpace 
      This provides the interfaces for discrete function spaces. 
      Discrete function spaces contain functions
      from a \ref FunctionSpaceInterface "function space"
      but the domain is defined by a grid
      or more precisly by a \ref GridPart "grid part".
 
      \remarks The interface for using a DiscreteFunctionSpace is
      defined by the class DiscreteFunctionSpaceInterface.

      @{
  */

  //! \brief enumerator for identification of spaces 
  enum DFSpaceIdentifier {  
    LagrangeSpace_id , //!< id for Lagrange Space 
    DGSpace_id ,       //!< id for Discontinuous Galerkin Space 
    CombinedSpace_id , //!< id for Combined Space 
    FiniteVolumeSpace_id , //!< id for Finite Volume Space 
    DFAdapter_id    //!< id for DiscreteFunctionSpace Adapter
  };



  //**************************************************************************
  //
  //  --DiscreteFunctionSpaceInterface
  //
  /**  
    \brief This is the interface for discrete function spaces. All methods
    declared here have to be implemented by the implementation class.

    The discrete function space always depends on a given grid. 
    For all diffrent element types of the grid the function space provides 
    a set of base functions for the different element types. 
    Because of the knowledge of on the one hand the grid an on the other
    hand the base functions sets, the discrete function space provides the size
    of the function space and a mapping from entity and local dof number
    to global dof number of the level of the entity.
    \note A DiscreteFunctionSpace is defined on a certain grid part.
   
    \interfaceclass 
  */
  template< class FunctionSpaceTraits >
  class DiscreteFunctionSpaceInterface
  : public FunctionSpaceTraits :: FunctionSpaceType
  {
  public:
    //- Typedefs and enums
    //! type of traits class 
    typedef FunctionSpaceTraits Traits; 
    //! type of DiscretefunctionSapce implementation (Barton-Nackman)
    typedef typename FunctionSpaceTraits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //! type of \ref Dune::FunctionSpaceInterface "function space" (define domain and range types)
    typedef typename FunctionSpaceTraits::FunctionSpaceType FunctionSpaceType;
    //! type of \ref Dune::BaseFunctionSetInterface "base function set" of this space 
    typedef typename FunctionSpaceTraits::BaseFunctionSetType BaseFunctionSetType;
    //! type of \ref Dune::DofMapperInterface "DoF mapper" of this space 
    typedef typename FunctionSpaceTraits::MapperType MapperType;
    //! type of underlying \ref GridPart "grid part" 
    typedef typename FunctionSpaceTraits::GridPartType GridPartType;
    //! type of underlying dune grid  
    typedef typename GridPartType::GridType GridType;
    //! type of used dune index set 
    typedef typename GridPartType::IndexSetType IndexSetType;
    
    /** \brief iterator type traversing the set of 
        entities defining the discrete  function space 
        (only codim 0 at the moment, to be revised)
    */
    typedef typename GridPartType:: template Codim<0>::IteratorType IteratorType;
    
  public:
    //- Public methods
    //! Constructor 
    DiscreteFunctionSpaceInterface ()
    : FunctionSpaceType() 
    {
    }

    //! Method provided by implementation

    /** \brief return type identifier of discrete function space 
        \return return type identifier of discrete function space
    */
    DFSpaceIdentifier type () const 
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().type());
      return asImp().type();
    }

    /** \brief get base function set for given entity
     *
     *  \param[in]  entity  entity for which base function is requested 
     *
     *  \returns BaseFunctionSet for the entity
     */
    template< class EntityType >
    inline const BaseFunctionSetType baseFunctionSet ( const EntityType &entity ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().baseFunctionSet( entity ) );
      return asImp().baseFunctionSet( entity );
    }
  
    /** \brief return true if the space contains globally continuous functions
     *
     *  For example, a LagrangeDiscreteFunctionSpace returns \b true and a
     *  DiscontinuousGalerkinSpace return \b false.
     *
     *  \return \b true if the space contians globally continous functions,
     *          \b false otherwise
     */
    inline bool continuous () const
    { 
      CHECK_INTERFACE_IMPLEMENTATION( asImp().continuous() );
      return asImp().continuous(); 
    }

    /** \brief get index of the sequence in grid sequences
     *
     *  \return number of current sequence
     */
    inline int sequence () const
    { 
      CHECK_INTERFACE_IMPLEMENTATION( asImp().sequence() );
      return asImp().sequence();
    }

    /** \brief get global order of space
     *
     *  \return order of space, i.e., the maximal polynomial order of base
     *          functions 
     */
    inline int order () const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION( asImp().order() );
      return asImp().order();
    } 
  
    /** \brief get a reference to the DoF mapper
     *
     *  \returns refernce to mapper
     */    
    inline MapperType &mapper () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().mapper() );
      return asImp().mapper();
    }
    
    /** \brief get reference to grid this discrete function space belongs to
     * 
     *  \hasdefault
     *
     *  \returns constant reference to grid  
     */
    inline const GridType &grid () const
    { 
      CHECK_INTERFACE_IMPLEMENTATION( asImp().grid() );
      return asImp().grid(); 
    }

    /** \brief get reference to grid this discrete function space belongs to
     * 
     *  \hasdefault
     *
     *  \returns reference to grid  
     */
    inline GridType &grid ()
    { 
      CHECK_INTERFACE_IMPLEMENTATION( asImp().grid() );
      return asImp().grid(); 
    }

    /** \brief get a reference to the associated grid partition
     *
     *  \hasdefault
     *
     *  \returns constant reference to the grid partition
     */
    inline const GridPartType &gridPart () const
    { 
      CHECK_INTERFACE_IMPLEMENTATION( asImp().gridPart() );
      return asImp().gridPart(); 
    }
    
    /** \brief get a reference to the associated grid partition
     *
     *  \hasdefault
     *
     *  \returns reference to the grid partition
     */
    inline GridPartType &gridPart ()
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().gridPart());
      return asImp().gridPart(); 
    }

    /** \brief Get a reference to the associated index set
     *
     *  \hasdefault
     *
     *  \returns const reference to index set
     */ 
    inline const IndexSetType &indexSet () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().indexSet() );
      return asImp().indexSet();
    }

    /** \brief get number of DoFs for this space
     *
     *  \hasdefault
     *
     *  \returns number of DoFs (degrees of freedom)
     */
    inline int size () const
    { 
      CHECK_INTERFACE_IMPLEMENTATION( asImp().size() );
      return asImp().size();
    }

    /** \brief get iterator pointing to the first entity of the associated grid
     *         partition
     *
     *  \hasdefault
     *
     *  \returns iterator pointing to first entity
     */
    inline IteratorType begin () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().begin() );
      return asImp().begin();
    }

    /** \brief get iterator pointing behind the last entity of the associated
     *         grid partition
     *
     *  \hasdefault
     *
     *  \returns iterator pointing behind last entity
     */
    inline IteratorType end () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().end() );
      return asImp().end();
    }

    /** \brief returns true if the grid has more than one geometry type
     *
     *  \hasdefault
     *
     *  \return \b true if the underlying grid has more than one geometry type
     *          (hybrid grid), \b false otherwise 
     */
    inline bool multipleGeometryTypes () const
    { 
      CHECK_INTERFACE_IMPLEMENTATION( asImp().multipleGeometryTypes() );
      return asImp().multipleGeometryTypes();
    }

    /** \brief returns true if base function sets depend on the entity
     *
     *  \hasdefault
     *
     *  \returns \b true if base function set depend on entities, \b false
     *           otherwise
     */
    inline bool multipleBaseFunctionSets () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().multipleBaseFunctionSets() );
      return asImp().multipleBaseFunctionSets();
    }
    
    /** \brief Map local DoF number to global DoF number
     *
     *  Maps an entity and a local DoF number to a global DoF number, i.e.,
     *  the index of the DoF within the DoF vector.
     *
     *  \hasdefault
     *
     *  \param[in]  entity    Entity for which mapping is done
     *  \param[in]  localDof  local dof number
     *
     *  \returns global DoF number
     */    
    template< class EntityType >
    inline int mapToGlobal ( const EntityType &entity,
                             const int localDof ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().mapToGlobal( entity, localDof ) );
      return asImp().mapToGlobal( entity, localDof );
    }

  protected:
    // Barton-Nackman trick 
    inline const DiscreteFunctionSpaceType &asImp () const
    { 
      return static_cast< const DiscreteFunctionSpaceType & >( *this );
    }

    // Barton-Nackman trick 
    inline DiscreteFunctionSpaceType &asImp()
    { 
      return static_cast< DiscreteFunctionSpaceType & >( *this ); 
    }
  }; // end class DiscreteFunctionSpaceInterface



  //---------------------------------------------------------------------------
  //-
  //-  --DiscreteFunctionSpaceDefault
  //-
  //-
  //---------------------------------------------------------------------------
  /** \brief This is the class with default implementations for discrete
     function.
     The  methods not marked with having a default
     in the interface class must be provided by
     the implementation; all other methods
     have a default implementation here.
     
     \remark An reference to the GridPart is
             stored in the default implementation. */
  template< class FunctionSpaceTraits >
  class DiscreteFunctionSpaceDefault
  : public DiscreteFunctionSpaceInterface< FunctionSpaceTraits >
  {
  public:
    typedef FunctionSpaceTraits Traits;
    
  private:
    typedef DiscreteFunctionSpaceDefault< Traits > ThisType;
    typedef DiscreteFunctionSpaceInterface< Traits > BaseType;
    
  public:
    typedef typename FunctionSpaceTraits :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    typedef typename FunctionSpaceTraits::GridPartType  GridPartType;
    typedef typename GridPartType:: template Codim<0>::IteratorType IteratorType;
    typedef typename GridPartType :: GridType GridType;
    typedef typename GridPartType :: IndexSetType IndexSetType;

    typedef TemporaryLocalFunctionFactory< DiscreteFunctionSpaceType >
      LocalFunctionFactoryType;
    typedef LocalFunctionStack< LocalFunctionFactoryType > LocalFunctionStorageType;

    typedef typename LocalFunctionStorageType :: LocalFunctionType
      LocalFunctionType;

  public:
    using BaseType :: mapper;
    
  protected:
    using BaseType :: asImp;

  protected:
    GridPartType &gridPart_;

    const LocalFunctionFactoryType lfFactory_;
    mutable LocalFunctionStorageType lfStorage_;
 
    // true if grid has more than one geometry type (hybrid grids)
    const bool multipleGeometryTypes_;
 
  public:
    //! constructor
    inline explicit DiscreteFunctionSpaceDefault( GridPartType &gridPart )
    : BaseType(),
      gridPart_( gridPart ),
      lfFactory_( asImp() ),
      lfStorage_( lfFactory_ ),
      multipleGeometryTypes_( AllGeomTypes< IndexSetType, GridType > 
                                :: multipleGeomTypes() )
    {
    }

    /** obtain a local function for an entity (to store intermediate values)
     *  
     *  \param[in]  entity  entity for which a local function is desired
     *
     *  \returns a local function backed by a small, fast array
     */
    template< class EntityType >
    LocalFunctionType localFunction ( const EntityType &entity ) const
    {
      return lfStorage_.localFunction( entity );
    }
    
    /** \copydoc Dune::DiscreteFunctionSpaceInterface::grid() const */
    inline const GridType &grid () const
    {
      return asImp().gridPart().grid();
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::grid() */
    inline GridType &grid ()
    {
      return asImp().gridPart().grid();
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::gridPart() const */
    inline const GridPartType &gridPart () const
    {
      return gridPart_;
    }
   
    /** \copydoc Dune::DiscreteFunctionSpaceInterface::gridPart() */
    inline GridPartType &gridPart ()
    {
      return gridPart_;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::indexSet() const */
    inline const IndexSetType &indexSet () const
    {
      return asImp().gridPart().indexSet();
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::size */
    inline int size () const
    {
      return mapper().size();
    }
  
    /** \copydoc Dune::DiscreteFunctionSpaceInterface::begin */
    inline IteratorType begin () const
    {
      return asImp().gridPart().template begin< 0 >();
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::end */
    inline IteratorType end () const
    {
      return asImp().gridPart().template end< 0 >();
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::multipleGeometryTypes */
    inline bool multipleGeometryTypes () const
    {
      return multipleGeometryTypes_;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::multipleBaseFunctionSets
     *
     *  \note The default implementation returns \b false.
     */
    inline bool multipleBaseFunctionSets () const
    {
      return false;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::mapToGlobal */
    template< class EntityType >
    inline int mapToGlobal ( const EntityType &entity,
                             const int localDof ) const
    {
      return mapper().mapToGlobal( entity, localDof );
    }
  };



  ////////////////////////////////////////////////////////////
  //
  //  DiscreteFunctionSpaceAdapter 
  //
  ////////////////////////////////////////////////////////////
  /** \brief Create Obejct that behaves like a discrete function space 
      without to provide functions with the iterator facilities.

      \note DiscreteFunctionSpaceAdapter is itself derived from the template
            argument FunctionSpaceImp. Hence, the constructor will call the
            default constructor of FunctionSpaceImp when this class is
            instanciated. So do not use discrete function spaces for the
            first template argument.
  */
  template< class FunctionSpaceImp, class GridPartImp >
  class DiscreteFunctionSpaceAdapter
  : public FunctionSpaceImp
  {
  public:
    // type of the underlying function space
    typedef FunctionSpaceImp FunctionSpaceType;
    //! type of the grid partition
    typedef GridPartImp GridPartType;

  private:
    typedef DiscreteFunctionSpaceAdapter< FunctionSpaceType, GridPartType >
      ThisType;
    typedef FunctionSpaceType BaseType;

  public:  
    enum { polynomialOrder = 111 };
   
    //! type of the grid
    typedef typename GridPartType :: GridType GridType;
    //! type of the index set 
    typedef typename GridPartType :: IndexSetType IndexSetType; 
    //! type of the grid iterator 
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType
      IteratorType;
    //- type of used entity
    typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
    
  protected:
    const GridPartType &gridPart_;

  public:
    //! constructor taking grid Part 
    inline explicit DiscreteFunctionSpaceAdapter ( const GridPartType &gridPart )
    : BaseType()
    , gridPart_( gridPart ) 
    {
    }

    //! copy constructor
    inline DiscreteFunctionSpaceAdapter( const ThisType &org )
    : gridPart_( org.gridPart_ ) 
    {
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::begin */
    inline IteratorType begin () const
    {
      return gridPart_.template begin< 0 >();
    }
    
    /** \copydoc Dune::DiscreteFunctionSpaceInterface::end */
    inline IteratorType end () const
    {
      return gridPart_.template end< 0 >();
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::gridPart */
    inline const GridPartType &gridPart () const
    {
      return gridPart_;
    }
    
    /** \copydoc Dune::DiscreteFunctionSpaceInterface::indexSet */
    inline const IndexSetType &indexSet () const
    {
      return gridPart_.indexSet();
    }
    
    /** \copydoc Dune::DiscreteFunctionSpaceInterface::grid */
    inline const GridType& grid () const
    {
      return gridPart_.grid();
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::continuous */
    inline bool continuous () const
    {
      return true;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::order */
    inline int order () const
    {
      return polynomialOrder;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::type */
    inline DFSpaceIdentifier type () const
    {
      return DFAdapter_id;
    }
  };

///@}  

  /**\ingroup HelperClasses 
     \brief 
     BaseFunctionSetSingletonFactory provides method createObject and
     deleteObject for the SingletonList  
  */
  template <class KeyImp, class ObjectImp, class ObjectFactoryImp>
  class BaseFunctionSetSingletonFactory
  { 
  public:
    //! create new BaseFunctionSet 
    static ObjectImp * createObject( const KeyImp & key )
    {
      ObjectFactoryImp fac(key); 
      return new ObjectImp(fac); 
    }
    
    //! delete BaseFunctionSet 
    static void deleteObject( ObjectImp * obj ) 
    {
      delete obj;
    }
  };
  
} // end namespace Dune 
#endif
