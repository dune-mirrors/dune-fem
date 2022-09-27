#ifndef DUNE_FEM_DISCRETEFUNCTIONSPACE_HH
#define DUNE_FEM_DISCRETEFUNCTIONSPACE_HH

// C++ includes
#include <cassert>

#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

// dune-common includes
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/dynvector.hh>

// dune-fem includes
#include <dune/fem/common/hybrid.hh>
#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/blockvectors/defaultblockvectors.hh>
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/misc/threads/threadsafevalue.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/auxiliarydofs.hh>
#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/version.hh>

// local includes
#include "allgeomtypes.hh"
#include "dofstorage.hh"


namespace Dune
{

  namespace Fem
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

    // ExportsDiscreteFunctionSpaceType
    // --------------------------------

    template< class T >
    class ExportsDiscreteFunctionSpaceType
    {
      typedef char Small;
      struct Big { char dummy[ 2 ]; };

      template< class U >
      static Small test ( const U &, typename U::DiscreteFunctionSpaceType * = nullptr );
      static Big test ( ... );

      static const T &makeT ();

      template< class U, bool >
      struct GetDiscreteFunctionSpaceType;

      template< class U >
      struct GetDiscreteFunctionSpaceType< U, true >
      {
        typedef typename U::DiscreteFunctionSpaceType Type;
      };

      template< class U >
      struct GetDiscreteFunctionSpaceType< U, false >
      {
        typedef void Type;
      };

    public:
      static const bool v = (sizeof( test( makeT() ) ) == sizeof( Small ));
      typedef typename GetDiscreteFunctionSpaceType< T, v >::Type Type;
    };



    // DFSpaceIdentifier
    // -----------------

    //! \brief enumerator for identification of spaces
    enum DFSpaceIdentifier {
      CombinedSpace_id,             //!< id for Combined Space
      DFAdapter_id,                 //!< id for DiscreteFunctionSpace Adapter
      DGSpace_id,                   //!< id for Discontinuous Galerkin Space
      FiniteVolumeSpace_id,         //!< id for Finite Volume Space
      FourierSpace_id,              //!< id for Fourier space
      GenericSpace_id,              //!< id for Generic Space
      LagrangeSpace_id,             //!< id for Lagrange Space
      RannacherTurekSpace_id,       //!< id for Rannacher-Turek space
      LegendreDGSpace_id,           //!< id for Legendre Discontinuous Galerkin Space
      HierarchicLegendreDGSpace_id, //!< id for Hierarchic Legendre Discontinuous Galerkin Space
      LagrangeDGSpace_id,           //!< id for Lagrange Discontinuous Galerkin Space
      LocalFiniteElementSpace_id         //!< id for local finite element space
    };

    inline std::string spaceName( const DFSpaceIdentifier id )
    {
      switch( id )
      {
        case CombinedSpace_id             : return "CombinedSpace";
        case DFAdapter_id                 : return "DFAdapter";
        case DGSpace_id                   : return "DiscontinuousGalerkinSpace";
        case FiniteVolumeSpace_id         : return "FiniteVolumeSpace";
        case FourierSpace_id              : return "FourierSpace";
        case GenericSpace_id              : return "GenericSpace";
        case LagrangeSpace_id             : return "LagrangeSpace";
        case RannacherTurekSpace_id       : return "RannacherTurekSpace";
        case LegendreDGSpace_id           : return "LegendreDGSpace";
        case HierarchicLegendreDGSpace_id : return "HierarchicLegendreDGSpace";
        case LagrangeDGSpace_id           : return "LagrangeDGSpace";
        default                           : return "unknown space";
      }
    }

    struct isGenericDiscreteFunctionSpace
    {};

    template< class DiscreteFunctionSpaceImp,
              class NewFunctionSpace>
    struct DifferentDiscreteFunctionSpace;

    template <class FunctionSpaceImp,
              class GridPartImp,
              int polOrd,
              class StorageImp,
              template <class,class,int,class> class DiscreteFunctionSpaceImp,
              class NewFunctionSpace>
    struct DifferentDiscreteFunctionSpace<
        DiscreteFunctionSpaceImp<FunctionSpaceImp,GridPartImp,polOrd,StorageImp>,
            NewFunctionSpace>
    {
      typedef DiscreteFunctionSpaceImp< NewFunctionSpace, GridPartImp, polOrd, StorageImp > Type;
    };

    template <class FunctionSpaceImp,
              class GridPartImp,
              int polOrd,
              bool caching,
              template <class,class,int,bool> class DiscreteFunctionSpaceImp,
              class NewFunctionSpace>
    struct DifferentDiscreteFunctionSpace<
        DiscreteFunctionSpaceImp<FunctionSpaceImp,GridPartImp,polOrd,caching>,
            NewFunctionSpace>
    {
      typedef DiscreteFunctionSpaceImp< NewFunctionSpace, GridPartImp, polOrd, caching > Type;
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
      //! type of \ref Dune::Fem::FunctionSpaceInterface "function space"
      typedef typename Traits :: FunctionSpaceType FunctionSpaceType;

    private:
      typedef FunctionSpaceType BaseType;

    public:
      //! type of \ref Dune::Fem::BasisFunctionSet "basis function set" of this space
      typedef typename Traits :: BasisFunctionSetType BasisFunctionSetType;
      //! type of block mapper of this space
      typedef typename Traits :: BlockMapperType BlockMapperType;

      typedef typename Traits::LocalBlockIndices LocalBlockIndices;

      //! size of local blocks
      static constexpr std::size_t localBlockSize = Hybrid::size( LocalBlockIndices() );

      //! type of underlying \ref GridPart "grid part"
      typedef typename Traits :: GridPartType GridPartType;

      //! type of underlying dune grid
      typedef typename GridPartType :: GridType GridType;
      //! type of used dune index set
      typedef typename GridPartType :: IndexSetType IndexSetType;
      /** \brief type of iterator for grid traversal

          \note Only grid traversal for codimension 0 is currently supported.
       */
      typedef typename GridPartType :: template Codim< Traits::codimension > :: IteratorType
        IteratorType;
      //! type of entity of codimension 0
      typedef typename GridPartType :: template Codim< Traits::codimension > :: EntityType EntityType;
      //! type of the intersections
      typedef typename GridPartType :: IntersectionType IntersectionType;
      //! type of auxiliary dofs
      typedef AuxiliaryDofs< GridPartType, BlockMapperType > AuxiliaryDofsType;

      //! deprecated type
      [[deprecated("Use AuxiliaryDofsType instead!")]]
      typedef AuxiliaryDofsType SlaveDofsType;

      /** \brief defines type of data handle for communication
          \param  DiscreteFunction  type of \ref Dune::Fem::DiscreteFunctionInterface
                                    "discrete function" to communicate
          \param  Operation         type of operation to perform on communication
                                    (defaults to copy)
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

      //! \brief typedef struct for defining the same discrete function space with a different function space
      template< class NewFunctionSpace >
      struct ToNewFunctionSpace
      {
        //! type of my discrete function space with new function space
        typedef typename DifferentDiscreteFunctionSpace< DiscreteFunctionSpaceType, NewFunctionSpace> :: Type Type;
      };

      //! \brief typedef struct for defining the same discrete function space with a different dimRange
      template< int newDimRange >
      struct ToNewDimRange
      {
        typedef typename ToNewDimRangeFunctionSpace< FunctionSpaceType, newDimRange > :: Type NewFunctionSpaceType;

        //! type of my discrete function space with new dim range
        typedef typename ToNewFunctionSpace< NewFunctionSpaceType > :: Type  Type;
      };

    private:
      static_assert( std::is_same<typename BaseType::DomainFieldType, typename GridType::ctype>::value,
                     "Domain field type of function space must equal field type of grid." );

    protected:
      DiscreteFunctionSpaceInterface ()
      {}

    public:

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

      /** \brief get basis function set for given entity

          \param[in]  entity  entity (of codim 0) for which base function is
                              requested

          \returns BasisFunctionSet for the entity
       */
      inline const BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().basisFunctionSet( entity ) );
        return asImp().basisFunctionSet( entity );
      }

      /** \brief returns true if the space contains only globally continuous
                 functions

          For example, a \ref Dune::Fem::LagrangeDiscreteFunctionSpace
          "Lagrange space" returns \b true while a \ref
          Dune::Fem::DiscontinuousGalerkinSpace "discontiuous Galerkin space" returns
          \b false.

          \returns \b true  if the space contians only globally continous
                            functions,
                   \b false otherwise
       */
      inline bool continuous () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().continuous() );
        return asImp().continuous();
      }

      /** \brief get index of the sequence in grid sequences

          \return number of current sequence
       */
      inline int sequence () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().sequence() );
        return asImp().sequence();
      }

      /** \brief get global order of space

          \return order of space, i.e., the maximal polynomial order of base
                  functions
       */
      inline int order () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().order() );
        return asImp().order();
      }

      /** \brief returns true if discrete functions over this space have zero jump
                 over the given intersection.

          For example, a \ref Dune::Fem::LagrangeDiscreteFunctionSpace
          "Lagrange space" returns \b true iff the intersection is conforming while a \ref
          Dune::Fem::DiscontinuousGalerkinSpace "discontiuous Galerkin space" always returns
          \b false.

          \param intersection Intersection for which we want to know the continuety
          \returns \b true  if the space contians functions which are continuous over the
                            intersection,
                   \b false otherwise
       */
      inline bool continuous (const IntersectionType &intersection) const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().continuous(intersection) );
        return asImp().continuous(intersection);
      }

      /** \brief get a reference to the block mapper

          \returns refernce to the block mapper
       */
      inline BlockMapperType &blockMapper () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().blockMapper() );
        return asImp().blockMapper();
      }

      /** \brief get reference to grid this discrete function space belongs to

          \returns constant reference to grid
       */
      inline const GridType &grid () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().grid() );
        return asImp().grid();
      }

      /** \brief get reference to grid this discrete function space belongs to

          \returns reference to grid
       */
      inline GridType &grid ()
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().grid() );
        return asImp().grid();
      }

      /** \brief get a reference to the associated grid partition

          \returns reference to the grid partition
       */
      inline GridPartType &gridPart ()
      {
        CHECK_INTERFACE_IMPLEMENTATION(asImp().gridPart());
        return asImp().gridPart();
      }

      /** \brief Get a reference to the associated index set

          \returns const reference to index set
       */
      inline const IndexSetType &indexSet () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().indexSet() );
        return asImp().indexSet();
      }

      /** \brief get number of DoFs for this space

          \returns number of DoFs (degrees of freedom)
       */
      inline int size () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().size() );
        return asImp().size();
      }

      /** \brief get number of primary DoFs for this space

          \returns number of primary DoFs (degrees of freedom that are owned by this process )
       */
      inline int primarySize () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().primarySize() );
        return asImp().primarySize();
      }

      /** \brief get number of auxiliary DoFs for this space

          \returns number of auxiliary DoFs (degrees of freedom that are NOT owned by this process)
       */
      inline int auxiliarySize () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().auxiliarySize() );
        return asImp().auxiliarySize();
      }

      /** \brief get iterator pointing to the first entity of the associated grid
                 partition

          \returns iterator pointing to first entity
       */
      inline IteratorType begin () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().begin() );
        return asImp().begin();
      }

      /** \brief get iterator pointing behind the last entity of the associated
                 grid partition

          \returns iterator pointing behind last entity
       */
      inline IteratorType end () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().end() );
        return asImp().end();
      }

      /** \brief apply a functor to each entity in the associated grid partition

          The functor must provide an the following operator
          \code
          template< class EntityType >
          void operator() ( const EntityType & );
          \endcode

          \param[in]  f  functor to apply
       */
      template< class FunctorType >
      inline void forEach ( FunctorType &f ) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().forEach( f ) );
      }

      /** \brief returns true if the grid has more than one geometry type

          \return \b true if the underlying grid has more than one geometry type
                  (hybrid grid), \b false otherwise
       */
      inline bool multipleGeometryTypes () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().multipleGeometryTypes() );
        return asImp().multipleGeometryTypes();
      }


      /** \brief returns true if base function sets depend on the entity

          \returns \b true if base function set depend on entities, \b false
                   otherwise
       */
      inline bool multipleBasisFunctionSets () const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().multipleBasisFunctionSets() );
        return asImp().multipleBasisFunctionSets();
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
                 default communication operation

          \param  discreteFunction  discrete function to be communicated
       */
      template <class DiscreteFunction>
      void communicate(DiscreteFunction& discreteFunction) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
            asImp().communicate( discreteFunction ) );
      }

      /** \brief communicate data for given discrete function

          \param      discreteFunction  discrete function to be communicated
          \param[in]  op                communication operation to use
                                        (see DFCommunicationOperation)
       */
      template <class DiscreteFunction, class Operation>
      void communicate(DiscreteFunction& discreteFunction, const Operation& op) const
      {
        CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
            asImp().communicate( discreteFunction , op ) );
      }

      /** \brief Creates DataHandle for given discrete function

          \param[in]  discreteFunction  \ref DiscreteFunctionInterface
                                        "discrete function" to create the data
                                        handle for
          \param[in]  operation         operation to perform on scatter
       */
      template< class DiscreteFunction, class Operation >
      inline typename CommDataHandle< DiscreteFunction, Operation > :: Type
      createDataHandle ( DiscreteFunction& discreteFunction,
                        const Operation &operation ) const
      {
        CHECK_INTERFACE_IMPLEMENTATION
          ( asImp().createDataHandle( discreteFunction, operation ) );
        return asImp().createDataHandle( discreteFunction, operation );
      }

      /** \brief get auxiliary dofs */
      const AuxiliaryDofsType& auxiliaryDofs() const
      {
        CHECK_INTERFACE_IMPLEMENTATION( asImp().auxiliaryDofs() );
        return asImp().auxiliaryDofs();
      }

      /** \brief deprecated method, use auxiliaryDofs */
      [[deprecated("Use auxiliaryDofs instead!")]]
      const AuxiliaryDofsType& slaveDofs() const
      {
        return auxiliaryDofs();
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
       \relates DiscreteFunctionSpaceInterface

        This is a default implemented equality operator for discrete function
        spaces. It assumes the mapper to be a singleton and then compares the
        addresses of the two mappers.

        Note that this method can be specialized by implementing another version
        that uses the exact traits of the discrete function space.
     */
    template< class Traits >
    inline bool operator== ( const DiscreteFunctionSpaceInterface< Traits > &X,
                             const DiscreteFunctionSpaceInterface< Traits > &Y )
    {
      return &(X.blockMapper()) == &(Y.blockMapper());
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
      : public DiscreteFunctionSpaceInterface< FunctionSpaceTraits >,
        public std::enable_shared_from_this< typename FunctionSpaceTraits::DiscreteFunctionSpaceType >
    {
    public:
      typedef FunctionSpaceTraits Traits;

    private:
      typedef DiscreteFunctionSpaceDefault< Traits > ThisType;
      typedef DiscreteFunctionSpaceInterface< Traits > BaseType;

    public:
      typedef typename Traits :: DiscreteFunctionSpaceType
        DiscreteFunctionSpaceType;

      typedef typename BaseType :: GridPartType GridPartType;
      typedef typename BaseType :: GridType GridType;
      typedef typename BaseType :: IndexSetType IndexSetType;
      typedef typename BaseType :: IteratorType IteratorType;
      typedef typename BaseType :: EntityType EntityType;

    protected:
      using BaseType :: asImp;

    public:
      using BaseType::localBlockSize;
      using BaseType :: blockMapper;
      using BaseType :: order ;

      //! type of DoF manager
      typedef DofManager< GridType > DofManagerType;

      //! type of communication manager
      typedef CommunicationManager< DiscreteFunctionSpaceType > CommunicationManagerType;

      typedef typename BaseType::BlockMapperType BlockMapperType;
      typedef typename BaseType::AuxiliaryDofsType AuxiliaryDofsType;

    protected:
      struct AuxiliaryDofsFactory
      {
        typedef std::pair< AuxiliaryDofsType, int > ObjectType;

        static ObjectType *createObject ( std::pair< GridPartType *, BlockMapperType * > key )
        {
          return new ObjectType( std::piecewise_construct, std::tie( *key.first, *key.second ), std::make_tuple( -1 ) );
        }

        static void deleteObject ( ObjectType *object ) { delete object; }
      };

      typedef SingletonList< std::pair< GridPartType *, BlockMapperType * >, std::pair< AuxiliaryDofsType, int >, AuxiliaryDofsFactory > AuxiliaryDofsProviderType;

    protected:
      GridPartType &gridPart_;

      typedef ThreadSafeValue< UninitializedObjectStack > LocalDofVectorStackType;
      typedef StackAllocator< typename BaseType::RangeFieldType, LocalDofVectorStackType* > LocalDofVectorAllocatorType;
      typedef Dune::DynamicVector< typename BaseType::RangeFieldType, LocalDofVectorAllocatorType > LocalDofVectorType;

      mutable LocalDofVectorStackType ldvStack_;
      mutable LocalDofVectorAllocatorType ldvAllocator_;

      typedef BasicTemporaryLocalFunction< ThisType, LocalDofVectorType > LocalFunctionType;

      // set of all geometry types possible
      typedef AllGeomTypes< IndexSetType, GridType > AllGeometryTypes;
      const AllGeometryTypes allGeomTypes_;

      // reference to dof manager
      DofManagerType& dofManager_;

      // communication manager
      const InterfaceType commInterface_;
      const CommunicationDirection commDirection_;
      mutable std::unique_ptr< CommunicationManagerType > communicator_;

    public:
      //! constructor
      explicit DiscreteFunctionSpaceDefault( GridPartType &gridPart,
          const InterfaceType commInterface = InteriorBorder_All_Interface,
          const CommunicationDirection commDirection = ForwardCommunication )
      : BaseType(),
        gridPart_( gridPart ),
        ldvStack_( 0 ),
        ldvAllocator_( &ldvStack_ ),
        allGeomTypes_( gridPart.indexSet() ),
        dofManager_( DofManagerType :: instance( gridPart.grid() ) ),
        commInterface_( commInterface ),
        commDirection_( commDirection )
      {}

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::sequence */
      inline int sequence () const
      {
        return dofManager_.sequence();
      }

      /** \brief default implementation of the method order
        *
        *  \return returns max polynomial order for each entity using the method order()
      */
      inline int order ( const EntityType& entity ) const
      {
        return asImp().basisFunctionSet( entity ).order();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::grid() const */
      inline const GridType &grid () const
      {
        return asImp().gridPart().grid();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::grid() */
      inline GridType &grid ()
      {
        return asImp().gridPart().grid();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::gridPart() const */
      inline GridPartType &gridPart () const
      {
        return gridPart_;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::indexSet() const */
      inline const IndexSetType &indexSet () const
      {
        return asImp().gridPart().indexSet();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::size */
      inline int size () const
      {
        return blockMapper().size() * localBlockSize ;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::primarySize */
      inline int primarySize () const
      {
        return auxiliaryDofs().primarySize() * localBlockSize;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::auxiliarySize */
      inline int auxiliarySize () const
      {
        // total size minus primary dofs
        return size() - primarySize();
      }

      //! \brief return the maximal number of dofs on entities
      inline int maxNumDofs () const
      {
        return blockMapper().maxNumDofs() * localBlockSize;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::begin() const
       *
       *  \note The default implementation uses the codim 0 iterators of the
       *        associated grid partition.
       */
      inline IteratorType begin () const
      {
        return asImp().gridPart().template begin< 0 >();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::end() const
       *
       *  \note The default implementation uses the codim 0 iterators of the
       *        associated grid partition.
       */
      inline IteratorType end () const
      {
        return asImp().gridPart().template end< 0 >();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::forEach(FunctorType &f) const

          \note The default implementation simply does the following:
          \code
          const IteratorType end = asImp().end();
          for( IteratorType it = asImp().begin(); it != end; ++it )
            f( *it );
          \endcode
       */
      template< class FunctorType >
      inline void forEach ( FunctorType &f ) const
      {
        const IteratorType end = asImp().end();
        for( IteratorType it = asImp().begin(); it != end; ++it )
          f( *it );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::multipleGeometryTypes */
      inline bool multipleGeometryTypes () const
      {
        return allGeomTypes_.multipleGeomTypes();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::multipleBasisFunctionSets
       *
       *  \note The default implementation returns \b false.
       */
      inline bool multipleBasisFunctionSets () const
      {
        return false;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::communicationInterface() */
      InterfaceType communicationInterface () const
      {
        return commInterface_;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::communicationInterface() */
      CommunicationDirection communicationDirection () const
      {
        return commDirection_;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::communicator() */
      const CommunicationManagerType& communicator() const
      {
        if( !communicator_ )
          communicator_.reset( new CommunicationManagerType( asImp(), commInterface_, commDirection_ ) );
        assert( communicator_ != 0 );
        return *communicator_;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::communicate(DiscreteFunction &discreteFunction) const */
      template <class DiscreteFunction>
      void communicate(DiscreteFunction& discreteFunction) const
      {
        // get type of default operation
        typedef typename DiscreteFunction :: DiscreteFunctionSpaceType :: template
          CommDataHandle< DiscreteFunction > :: OperationType  DefaultOperationType;

        // default communication operation
        DefaultOperationType operation;

        // exchange data
        communicate( discreteFunction, operation );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::communicate(DiscreteFunction &discreteFunction, const Operation &operation) const */
      template <class DiscreteFunction, class Operation>
      void communicate(DiscreteFunction& discreteFunction, const Operation& op ) const
      {
        if constexpr ( std::is_base_of< IsDiscreteFunction, DiscreteFunction >::value )
        {
          static_assert( std::is_same< typename DiscreteFunctionSpaceType::BlockMapperType,
              typename DiscreteFunction::DiscreteFunctionSpaceType::BlockMapperType >::value &&
              localBlockSize == static_cast< std::size_t >(DiscreteFunction::DiscreteFunctionSpaceType::localBlockSize),
              "DiscreteFunctionSpaceDefault::communicate cannot be called with discrete functions defined over a different space" );
          communicator().exchange( discreteFunction, op );
        }
        else
        {
          static_assert( std::is_base_of< IsBlockVector, DiscreteFunction> :: value, "DiscreteFunctionSpaceDefault::communicate needs at least a BlockVectorInterface and derived");
          communicator().exchange( asImp(), discreteFunction, op );
        }
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::createDataHandle(DiscreteFunction &discreteFunction.const Operation &operation) const

          \note The default implementation is
          \code
          return CommDataHandle< DiscreteFunction, Operation > :: Type( discreteFunction );
          \endcode
       */
      template< class DiscreteFunction, class Operation >
      inline typename BaseType
        :: template CommDataHandle< DiscreteFunction, Operation > :: Type
      createDataHandle( DiscreteFunction &discreteFunction,
                        const Operation& operation ) const
      {
        static_assert( std::is_same< typename DiscreteFunctionSpaceType::BlockMapperType,
            typename DiscreteFunction::DiscreteFunctionSpaceType::BlockMapperType >::value &&
            localBlockSize == static_cast< std::size_t >(DiscreteFunction::DiscreteFunctionSpaceType::localBlockSize),
            "DiscreteFunctionSpaceDefault::createDataHandle cannot be called with discrete functions defined over a different space" );
        return typename BaseType
          :: template CommDataHandle< DiscreteFunction, Operation >
          :: Type( discreteFunction, operation );
      }

      /** \brief get auxiliary dofs */
      const AuxiliaryDofsType& auxiliaryDofs() const
      {
        if ( !auxiliaryDofs_ )
          auxiliaryDofs_.reset( &(AuxiliaryDofsProviderType::getObject( std::make_pair( &this->gridPart(), &this->blockMapper() ) )) );
        const int sequence = asImp().sequence();
        if( auxiliaryDofs_->second != sequence )
        {
          auxiliaryDofs_->first.rebuild();
          auxiliaryDofs_->second = sequence;
        }
        return auxiliaryDofs_->first;
      }

      /** \brief default implementation of addFunction does nothing at the moment */
      template <class DiscreteFunction>
      void addFunction( DiscreteFunction& df ) const
      {
      }

      /** \brief default implementation of removeFunction does nothing at the moment */
      template <class DiscreteFunction>
      void removeFunction( DiscreteFunction& df ) const
      {
      }

      /** \brief default implementation of adapt does nothing,
                 its only used in PAdaptiveLagrangeSpace */
      template <class Vector>
      void adapt( const Vector& polynomialOrders, const int polOrderShift = 0 ) const
      {
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
      mutable std::unique_ptr< std::pair< AuxiliaryDofsType, int >, typename AuxiliaryDofsProviderType::Deleter > auxiliaryDofs_;
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
    , public std::enable_shared_from_this< DiscreteFunctionSpaceAdapter<FunctionSpaceImp,GridPartImp> >
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
      enum { polynomialOrder = 6 }; // default polynomial order basically for determination of quadrature orders

      //! type of the grid
      typedef typename GridPartType :: GridType GridType;
      //! type of the index set
      typedef typename GridPartType :: IndexSetType IndexSetType;
      //! type of the grid iterator
      typedef typename GridPartType :: template Codim< 0 > :: IteratorType
        IteratorType;
      //- type of used entity
      typedef typename GridType :: template Codim< 0 > :: Entity EntityType;
      //- type of intersections
      typedef typename GridPartType :: IntersectionType IntersectionType;

      //! type of communication manager (only the default communication is valid here)
      typedef DefaultCommunicationManager< ThisType > CommunicationManagerType;

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

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::begin */
      inline IteratorType begin () const
      {
        return gridPart_.template begin< 0 >();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::end */
      inline IteratorType end () const
      {
        return gridPart_.template end< 0 >();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::forEach(FunctorType &f) const */
      template< class FunctorType >
      inline void forEach ( FunctorType &f ) const
      {
        const IteratorType endit = end();
        for( IteratorType it = begin(); it != endit; ++it )
          f( *it );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::gridPart */
      inline const GridPartType &gridPart () const
      {
        return gridPart_;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::indexSet */
      inline const IndexSetType &indexSet () const
      {
        return gridPart_.indexSet();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::grid */
      inline const GridType& grid () const
      {
        return gridPart_.grid();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      inline bool continuous () const
      {
        return true;
      }
      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      inline bool continuous (const IntersectionType &intersection) const
      {
        return true;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      inline int order () const
      {
        return order_;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      inline int order ( const EntityType& ) const
      {
        return order();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      inline DFSpaceIdentifier type () const
      {
        return DFAdapter_id;
      }
    };

  ///@}

    /**\ingroup HelperClasses
       \brief
       BasisFunctionSetSingletonFactory provides method createObject and
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

  } // namespace Fem

} // namespace Dune
#endif // #ifndef DUNE_FEM_DISCRETEFUNCTIONSPACE_HH
