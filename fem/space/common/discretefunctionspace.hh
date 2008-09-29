#ifndef DUNE_DISCRETEFUNCTIONSPACE_HH
#define DUNE_DISCRETEFUNCTIONSPACE_HH

//- system includes
#include <assert.h>

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/fem/space/common/geometryconversion.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/basefunctioninterface.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/function/localfunction/localfunctionwrapper.hh>
#include <dune/fem/function/localfunction/temporarylocalfunction.hh>
#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/communicationmanager.hh>


//- local includes 
#include "allgeomtypes.hh"
#include "dofstorage.hh"

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
    //! type of traits class 
    typedef FunctionSpaceTraits Traits;

    //! type of DiscretefunctionSapce implementation (Barton-Nackman)
    typedef typename Traits :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    //! type of \ref Dune::FunctionSpaceInterface "function space"
    typedef typename Traits :: FunctionSpaceType FunctionSpaceType;
    
  private:
    typedef FunctionSpaceType BaseType;

  public:
    //! type of \ref Dune::BaseFunctionSetInterface "base function set" of this space 
    typedef typename Traits :: BaseFunctionSetType BaseFunctionSetType;
    //! type of \ref Dune::DofMapper "DoF mapper" of this space
    typedef typename Traits :: MapperType MapperType;
    //! type of block mapper of this space
    typedef typename Traits :: BlockMapperType BlockMapperType;

    //! type of underlying \ref GridPart "grid part" 
    typedef typename Traits :: GridPartType GridPartType;

    //! type of underlying dune grid  
    typedef typename GridPartType :: GridType GridType;
    //! type of used dune index set 
    typedef typename GridPartType :: IndexSetType IndexSetType;
    /** \brief type of iterator for grid traversal
     *
     *  \note Only grid traversal for codimension 0 is currently supported.
     */
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType
      IteratorType;
    //! type of entity of codimension 0
    typedef typename IteratorType :: Entity EntityType;
 
    /** \brief defines type of data handle for communication
     *  \param  DiscreteFunction  type of \ref Dune::DiscreteFunctionInterface
     *                            "discrete function" to communicate
     *  \param  Operation         type of operation to perform on communication
     *                            (defaults to copy)
     */
    template< class DiscreteFunction,
              class Operation =  // get default type from traits 
                typename Traits :: template CommDataHandle< DiscreteFunction > :: OperationType
            >
    struct CommDataHandle
    {
      //! type of communication data handle
      typedef typename Traits
        :: template CommDataHandle< DiscreteFunction, Operation > :: Type
        Type;

      //! type of operation to perform on scatter
      typedef typename Traits
        :: template CommDataHandle< DiscreteFunction, Operation > :: OperationType
        OperationType;
    };

    //! type of communication manager 
    typedef CommunicationManager< DiscreteFunctionSpaceType > CommunicationManagerType;

  private:
    dune_static_assert( (Conversion<typename BaseType::DomainFieldType,
                                    typename GridType::ctype>::sameType),
                        "Domain field type of function space must equal field type of grid." );
  public:
    //! Constructor 
    inline DiscreteFunctionSpaceInterface ()
    : BaseType()
    {}

    // Methods Provided by the Implementation
    // --------------------------------------

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
     *  \param[in]  entity  entity (of codim 0) for which base function is
     *                      requested
     *
     *  \returns BaseFunctionSet for the entity
     */
    inline const BaseFunctionSetType baseFunctionSet ( const EntityType &entity ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().baseFunctionSet( entity ) );
      return asImp().baseFunctionSet( entity );
    }
  
    /** \brief returns true if the space contains DoFs of given codimension
     * 
     *  \param[in]  codim  codimension to check for DoFs
     *  
     *  \returns \b true if codimension contains DoFs,
     *           \b false otherwise
     */
    inline bool contains ( const int codim ) const
    { 
      CHECK_INTERFACE_IMPLEMENTATION( asImp().contains( codim ) );
      return asImp().contains( codim ); 
    }

    /** \brief returns true if the space contains only globally continuous
     *         functions
     *
     *  For example, a \ref Dune::LagrangeDiscreteFunctionSpace
     *  "Lagrange space" returns \b true while a \ref
     *  Dune::DiscontinuousGalerkinSpace "discontiuous Galerkin space" returns
     *  \b false.
     *
     *  \returns \b true  if the space contians only globally continous
     *                    functions,
     *           \b false otherwise
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

    /** \brief get a reference to the block mapper
     *
     *  \returns refernce to the block mapper
     */    
    inline BlockMapperType &blockMapper () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().blockMapper() );
      return asImp().blockMapper();
    }
    
    /** \brief get reference to grid this discrete function space belongs to
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
     *  \returns reference to grid  
     */
    inline GridType &grid ()
    { 
      CHECK_INTERFACE_IMPLEMENTATION( asImp().grid() );
      return asImp().grid(); 
    }

    /** \brief get a reference to the associated grid partition
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
     *  \returns reference to the grid partition
     */
    inline GridPartType &gridPart ()
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().gridPart());
      return asImp().gridPart(); 
    }

    /** \brief Get a reference to the associated index set
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
     *  \returns iterator pointing behind last entity
     */
    inline IteratorType end () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().end() );
      return asImp().end();
    }
    
    /** \brief apply a functor to each entity in the associated grid partition
     *
     *  The functor must provide an the following operator
     *  \code
     *  template< class EntityType >
     *  void operator() ( const EntityType & );
     *  \endcode
     *
     *  \param[in]  f  functor to apply
     */
    template< class FunctorType >
    inline void forEach ( FunctorType &f ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().forEach( f ) );
    }

    /** \brief returns true if the grid has more than one geometry type
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
     *  \param[in]  entity    entity (of codim 0) for which the mapping is done
     *  \param[in]  localDof  local dof number
     *
     *  \returns global DoF number
     */    
    inline int mapToGlobal ( const EntityType &entity,
                             const int localDof ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().mapToGlobal( entity, localDof ) );
      return asImp().mapToGlobal( entity, localDof );
    }


    /** \brief return the communication interface appropriate for this space 
        \return communication interface 
    */
    InterfaceType communicationInterface() const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().communicationInterface() );
      return asImp().communicationInterface();
    }

    /** \brief return the communication direction appropriate for this space 
        \return communication direction  
    */
    CommunicationDirection communicationDirection() const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().communicationDirection() );
      return asImp().communicationDirection();
    }

    /** \brief return reference to communicator (see CommunicationManager)
        \return reference to communicator 
    */
    const CommunicationManagerType& communicator() const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().communicator() );
      return asImp().communicator();
    }

    /** \brief communicate data for given discrete function using the space's
     *         default communication operation
     *
     *  \param  discreteFunction  discrete function to be communicated
     */
    template <class DiscreteFunction>
    void communicate(DiscreteFunction& discreteFunction) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
          asImp().communicate( discreteFunction ) );
    }

    /** \brief communicate data for given discrete function
     *
     *  \param      discreteFunction  discrete function to be communicated
     *  \param[in]  op                communication operation to use
     *                                (see DFCommunicationOperation)
     */
    template <class DiscreteFunction, class Operation>
    void communicate(DiscreteFunction& discreteFunction, const Operation* op) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
          asImp().communicate( discreteFunction , op ) );
    }

    /** \brief return maximal number of local DoFs
     *
     *  \returns an upper bound for the number of local DoFs
     */
    inline int maxNumLocalDofs () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().maxNumDofs() );
      return asImp().maxNumDofs();
    }

    /** \brief Creates DataHandle for given discrete function
     *
     *  \param[in]  discreteFunction  \ref DiscreteFunctionInterface
     *                                "discrete function" to create the data
     *                                handle for
     *  \param[in]  operation         operation to perform on scatter
     */
    template< class DiscreteFunction, class Operation >
    inline typename CommDataHandle< DiscreteFunction, Operation > :: Type
    createDataHandle ( DiscreteFunction& discreteFunction,
                      const Operation *operation ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION
        ( asImp().createDataHandle( discreteFunction, operation ) );
      return asImp().createDataHandle( discreteFunction, operation );
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



  /** \brief check two spaces for equality
   *  \relates DiscreteFunctionSpaceInterface
   *  
   *  This is a default implemented equality operator for discrete function
   *  spaces. It assumes the mapper to be a singleton and then compares the
   *  addresses of the two mappers.
   *
   *  Note that this method can be specialized by implementing another version
   *  that uses the exact traits of the discrete function space.
   */
  template< class Traits >
  inline bool operator== ( const DiscreteFunctionSpaceInterface< Traits > &X,
                           const DiscreteFunctionSpaceInterface< Traits > &Y )
  {
    return &(X.mapper()) == &(Y.mapper());
  }



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

    typedef typename BaseType :: GridPartType GridPartType;
    typedef typename BaseType :: GridType GridType;
    typedef typename BaseType :: IndexSetType IndexSetType;
    typedef typename BaseType :: IteratorType IteratorType;
    typedef typename BaseType :: EntityType EntityType;

    typedef TemporaryLocalFunctionFactory< DiscreteFunctionSpaceType >
      LocalFunctionFactoryType;
    typedef LocalFunctionStack< LocalFunctionFactoryType > LocalFunctionStorageType;

    typedef typename LocalFunctionStorageType :: LocalFunctionType
      LocalFunctionType;

  protected:
    using BaseType :: asImp;

  public:
    using BaseType :: mapper;

    //! type of DoF manager
    typedef DofManager< GridType > DofManagerType;
    //! type of DoF manager factory
    typedef DofManagerFactory< DofManagerType > DofManagerFactoryType;

    //! type of communication manager 
    typedef CommunicationManager< DiscreteFunctionSpaceType > CommunicationManagerType;
  protected:
    mutable GridPartType &gridPart_;

    const LocalFunctionFactoryType lfFactory_;
    mutable LocalFunctionStorageType lfStorage_;
 
    // set of all geometry types possible 
    typedef AllGeomTypes< IndexSetType, GridType > AllGeometryTypes;
    const AllGeometryTypes allGeomTypes_;

    // reference to dof manager 
    DofManagerType& dofManager_;

    // communication manager 
    const InterfaceType commInterface_;
    const CommunicationDirection commDirection_;
    mutable CommunicationManagerType *communicator_;
 
  public:
    //! constructor
    explicit DiscreteFunctionSpaceDefault( GridPartType &gridPart,
        const InterfaceType commInterface = InteriorBorder_All_Interface,
        const CommunicationDirection commDirection = ForwardCommunication )
    : BaseType(),
      gridPart_( gridPart ),
      lfFactory_( asImp() ),
      lfStorage_( lfFactory_ ),
      allGeomTypes_( gridPart.indexSet() ),
      dofManager_( DofManagerFactoryType :: getDofManager( gridPart.grid() ) ),
      commInterface_( commInterface ),
      commDirection_( commDirection ),
      communicator_( 0 )
    {}

  protected:
    ~DiscreteFunctionSpaceDefault ()
    {
      if( communicator_ != 0 )
        delete communicator_;
    }

  public:
    /** \copydoc Dune::DiscreteFunctionSpaceInterface::sequence */
    inline int sequence () const
    { 
      return dofManager_.sequence();
    }

    /** obtain a local function for an entity (to store intermediate values)
     *  
     *  \param[in]  entity  entity (of codim 0) for which a local function is
     *                      desired
     *
     *  \returns a local function backed by a small, fast array
     */
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
    inline GridPartType &gridPart () const
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
  
    /** \copydoc Dune::DiscreteFunctionSpaceInterface::begin() const
     *
     *  \note The default implementation uses the codim 0 iterators of the
     *        associated grid partition.
     */
    inline IteratorType begin () const
    {
      return asImp().gridPart().template begin< 0 >();
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::end() const
     *
     *  \note The default implementation uses the codim 0 iterators of the
     *        associated grid partition.
     */
    inline IteratorType end () const
    {
      return asImp().gridPart().template end< 0 >();
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::forEach(FunctorType &f) const
     *
     *  \note The default implementation simply does the following:
     *  \code
     *  const IteratorType end = asImp().end();
     *  for( IteratorType it = asImp().begin(); it != end; ++it )
     *    f( *it );
     *  \endcode
     */
    template< class FunctorType >
    inline void forEach ( FunctorType &f ) const
    {
      const IteratorType end = asImp().end();
      for( IteratorType it = asImp().begin(); it != end; ++it )
        f( *it );
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::multipleGeometryTypes */
    inline bool multipleGeometryTypes () const
    {
      return allGeomTypes_.multipleGeomTypes();
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::multipleBaseFunctionSets
     *
     *  \note The default implementation returns \b false.
     */
    inline bool multipleBaseFunctionSets () const
    {
      return false;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::mapToGlobal(const EntityType &entity,const int localDof) const */
    inline int mapToGlobal ( const EntityType &entity,
                             const int localDof ) const
    {
      return mapper().mapToGlobal( entity, localDof );
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::maxNumLocalDofs */
    inline int maxNumLocalDofs () const
    {
      return mapper().maxNumDofs();
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::communicationInterface() */
    InterfaceType communicationInterface () const
    {
      return commInterface_;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::communicationInterface() */
    CommunicationDirection communicationDirection () const
    {
      return commDirection_;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::communicator() */
    const CommunicationManagerType& communicator() const
    {
      if( communicator_ == 0 )
      {
        communicator_
          = new CommunicationManagerType( asImp(), commInterface_, commDirection_ );
      }
      assert( communicator_ != 0 );
      return *communicator_;
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::communicate(DiscreteFunction &discreteFunction) const */
    template <class DiscreteFunction>
    void communicate(DiscreteFunction& discreteFunction) const
    {
      // get type of default operation 
      typedef typename DiscreteFunction :: DiscreteFunctionSpaceType :: template
        CommDataHandle< DiscreteFunction > :: OperationType  DefaultOperationType;

      // exchange data 
      communicate( discreteFunction, (DefaultOperationType*) 0);
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::communicate(DiscreteFunction &discreteFunction, const Operation *) const */
    template <class DiscreteFunction, class Operation>
    void communicate(DiscreteFunction& discreteFunction, const Operation *op ) const
    {
      communicator().exchange( discreteFunction, (Operation *) 0);
    }

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::createDataHandle(DiscreteFunction &discreteFunction.const Operation *operation) const
     *
     *  \note The default implementation is
     *  \code
     *  return CommDataHandle< DiscreteFunction, Operation > :: Type( discreteFunction );
     *  \endcode
     */
    template< class DiscreteFunction, class Operation >
    inline typename BaseType
      :: template CommDataHandle< DiscreteFunction, Operation > :: Type
    createDataHandle( DiscreteFunction &discreteFunction,
                      const Operation *operation ) const
    {
      return typename BaseType
        :: template CommDataHandle< DiscreteFunction, Operation >
        :: Type( discreteFunction );
    }

  protected:  
    /** \brief returns true if the grid has more than one geometry type
     *
     *  \return \b true if the underlying grid has more than one geometry type
     *          (hybrid grid), \b false otherwise 
     */
    inline const std::vector<GeometryType>& geomTypes(int codim) const 
    { 
      return allGeomTypes_.geomTypes(codim);
    }

    // only combined space should use geomTypes 
    template <class , int , DofStoragePolicy> friend class CombinedSpace;
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
    const unsigned int order_;

  public:
    //! constructor taking grid Part 
    inline explicit DiscreteFunctionSpaceAdapter
      ( const GridPartType &gridPart,
        unsigned int order = polynomialOrder )
    : BaseType(),
      gridPart_( gridPart ),
      order_( order )
    {
    }

    //! copy constructor
    inline DiscreteFunctionSpaceAdapter( const ThisType &other )
    : BaseType( other ),
      gridPart_( other.gridPart_ ),
      order_( other.order_ )
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
      return order_;
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
