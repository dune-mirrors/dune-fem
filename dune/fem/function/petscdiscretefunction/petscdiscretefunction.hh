#ifndef DUNE_FEM_PETSCDISCRETEFUNCTION_HH
#define DUNE_FEM_PETSCDISCRETEFUNCTION_HH

#include <memory>
#include <string>
#include <utility>

#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/vectorfunction.hh>

#if HAVE_PETSC

#include <dune/common/dynvector.hh>
#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/functor.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/misc/petsc/petscvector.hh>

#include <dune/fem/space/common/restrictprolonginterface.hh>

namespace Dune
{
  namespace Fem
  {

    template< class DiscreteFunctionSpace >
    class PetscDiscreteFunction;

    template< class DofProxy, class Allocator >
    struct AssignVectorReference< Dune::DynamicVector< DofProxy, Allocator >  >
    {
      // we need to overload this
      typedef Dune::DynamicVector< DofProxy, Allocator > Vector;

      AssignVectorReference ( Vector &vector )
        : vector_( vector )
      {}

      void operator() ( const std::size_t index, DofProxy value ) const
      {
        vector_[ index ].assign( value );
      }

    protected:
      Vector &vector_;
    };



    /** \class DiscreteFunctionTraits for PetscDiscreteFunction
     *  \brief Traits class for a DiscreteFunction
     *
     *  \tparam  DiscreteFunctionSpace   space the discrete function lives in
     *  \tparam  DofVector             implementation class of the block vector
     */
    template< typename DiscreteFunctionSpace >
    struct DiscreteFunctionTraits< PetscDiscreteFunction< DiscreteFunctionSpace > >
      : public DefaultDiscreteFunctionTraits< DiscreteFunctionSpace, PetscVector< DiscreteFunctionSpace > >
    {
      typedef PetscVector< DiscreteFunctionSpace >  DofVectorType;
      typedef PetscDiscreteFunction< DiscreteFunctionSpace > DiscreteFunctionType;

      typedef typename DofVectorType::DofBlockType          DofBlockType;
      typedef typename DofBlockType::DofProxy DofProxyType;

      typedef ThreadSafeValue< UninitializedObjectStack > LocalDofVectorStackType;
      typedef StackAllocator< DofProxyType, LocalDofVectorStackType* > LocalDofVectorAllocatorType;
      typedef Dune::DynamicVector< DofProxyType, LocalDofVectorAllocatorType > LocalDofVectorType;

      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };



    // PetscDiscreteFunction
    // ---------------------

    template <class DiscreteFunctionSpace>
    class PetscDiscreteFunction
    : public DiscreteFunctionDefault< PetscDiscreteFunction< DiscreteFunctionSpace > >
    {
      typedef PetscDiscreteFunction< DiscreteFunctionSpace > ThisType;
      typedef DiscreteFunctionDefault< ThisType > BaseType;

    public:
      typedef typename BaseType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
      typedef typename BaseType :: DofVectorType              DofVectorType;

      // generic assign method
      using BaseType::assign;

      /** \brief Constructor to use if the vector storing the dofs does not exist yet
       *
       *  \param[in]  name         name of the discrete function
       *  \param[in]  space        space the discrete function lives in
       */
      PetscDiscreteFunction( const std::string& name,
                             const DiscreteFunctionSpaceType& space )
        : BaseType( name, space ),
          memObject_(),
          dofVector_( allocateDofStorage( space ) )
      {}

      /** \brief Constructor to use if the vector storing the dofs already exists
       *
       *  \param[in]  name         name of the discrete function
       *  \param[in]  space        space the discrete function lives in
       *  \param[in]  dofVector    reference to the dof vector
       */
      PetscDiscreteFunction( const std::string& name,
                             const DiscreteFunctionSpaceType& space,
                             DofVectorType& dofVector )
        : BaseType( name, space ),
          memObject_(),
          dofVector_( dofVector )
      {}

      /** \brief Copy constructor */
      PetscDiscreteFunction( const ThisType& other )
        : BaseType( "copy of " + other.name(), other.space() ),
          memObject_(),
          dofVector_( allocateDofStorage( other.space() ) )
      {
        assign( other );
      }

      /** \brief Move constructor */
      PetscDiscreteFunction( ThisType&& other )
        : BaseType( static_cast< BaseType && >( other ) ),
          memObject_( std::move( other.memObject_ ) ),
          dofVector_( other.dofVector_ )
      {}

      PetscDiscreteFunction () = delete;
      ThisType& operator= ( const ThisType& ) = delete;
      ThisType& operator= ( ThisType&& ) = delete;

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::enableDofCompression() */
      void enableDofCompression ()
      {
        if( memObject_ )
          memObject_->enableDofCompression();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::assign(const DiscreteFunctionInterfaceType &g)
       *  \note This is a specialization when the right hand side is an AdaptiveDiscreteFunction */
      void assign( const AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > &g )
      {
        // call more efficient assign on PetscVector
        dofVector().assignVector( g.dofVector() );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::assign(const DiscreteFunctionInterfaceType &g)
       *  \note This is a specialization when the right hand side is an ISTLBlockVectorDiscreteFunction */
      template < class Block >
      void assign( const ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType, Block > &g )
      {
        // call more efficient assign on PetscVector
        dofVector().assignVector( g.dofVector() );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::assign(const DiscreteFunctionInterfaceType &g)
       *  \note This is a specialization when the right hand side is a VectorDiscreteFunction */
      template < class Vector >
      void assign( const VectorDiscreteFunction< DiscreteFunctionSpaceType, Vector > &g )
      {
        // call more efficient assign on PetscVector
        dofVector().assignVector( g.dofVector() );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::beginAssemble() */
      template< class AssembleOperation >
      void beginAssemble ()
      {
        BaseType :: template beginAssemble< AssembleOperation > ();
        dofVector().beginAssemble();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::endAssemble() */
      template< class AssembleOperation >
      void endAssemble ()
      {
        BaseType :: template endAssemble< AssembleOperation > ();
        dofVector().endAssemble();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dofVector() */
      DofVectorType& dofVector() { return dofVector_; }
      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dofVector() */
      const DofVectorType& dofVector() const { return dofVector_; }

      /** \brief obtain a constant pointer to the underlying PETSc Vec */
      const Vec* petscVec () const { return dofVector().getVector(); }

      /** \brief obtain a pointer to the underlying PETSc Vec */
      Vec* petscVec () { return dofVector().getVector(); }

    protected:
      typedef typename DiscreteFunctionSpaceType :: BlockMapperType  BlockMapperType;
      typedef PetscManagedDofStorage< DiscreteFunctionSpaceType, BlockMapperType > PetscManagedDofStorageType;

      // allocate managed dof storage
      DofVectorType& allocateDofStorage ( const DiscreteFunctionSpaceType &space )
      {
        memObject_.reset( new PetscManagedDofStorageType( space, space.blockMapper() ) );

        return memObject_->getArray();
      }

      // pointer to allocated DofVector
      std::unique_ptr< PetscManagedDofStorageType > memObject_;

      // dof vector impl
      DofVectorType& dofVector_;
    };


    template <class DiscreteFunctionSpace>
    class RestrictProlongDefault< PetscDiscreteFunction< DiscreteFunctionSpace > >
      : public RestrictProlongInterfaceDefault< RestrictProlongTraits< RestrictProlongDefault< PetscDiscreteFunction< DiscreteFunctionSpace > >,
                                                typename DiscreteFunctionSpace::DomainFieldType > >
    {
      typedef PetscDiscreteFunction< DiscreteFunctionSpace >   DiscreteFunction;
      typedef RestrictProlongDefault< DiscreteFunction > ThisType;
      typedef RestrictProlongInterfaceDefault< RestrictProlongTraits< ThisType, typename DiscreteFunction::DomainFieldType > > BaseType;

      typedef AdaptiveDiscreteFunction< DiscreteFunctionSpace >         AdaptiveDiscreteFunctionType;
      typedef RestrictProlongDefault < AdaptiveDiscreteFunctionType >   AdaptiveRestrictProlongType;

    public:
      typedef DiscreteFunction DiscreteFunctionType;

      typedef typename BaseType::DomainFieldType DomainFieldType;

      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
      typedef typename DiscreteFunctionType::GridPartType GridPartType;

      typedef DefaultLocalRestrictProlong< DiscreteFunctionSpaceType > LocalRestrictProlongType;

      explicit RestrictProlongDefault ( DiscreteFunctionType &discreteFunction )
      : discreteFunction_( discreteFunction ),
        adaptiveFunction_( discreteFunction_.name()+"-adaptive", discreteFunction_.space() ),
        rpOp_( adaptiveFunction_ ),
        initialized_( false )
      {
      }

      RestrictProlongDefault ( const RestrictProlongDefault& other )
      : discreteFunction_( const_cast< DiscreteFunction& > (other.discreteFunction_) ),
        adaptiveFunction_( other.adaptiveFunction_.name()+"-copy", discreteFunction_.space() ),
        rpOp_( adaptiveFunction_ ),
        initialized_( false )
      {
      }

      void setFatherChildWeight ( const DomainFieldType &weight ) const
      {
        rpOp_.setFatherChildWeight( weight );
      }

      //! restrict data to father
      template< class Entity >
      void restrictLocal ( const Entity &father, const Entity &son, bool initialize ) const
      {
        assert( initialized_ );
        rpOp_.restrictLocal( father, son, initialize );
      }

      //! prolong data to children
      template< class Entity >
      void prolongLocal ( const Entity &father, const Entity &son, bool initialize ) const
      {
        assert( initialized_ );
        rpOp_.prolongLocal( father, son, initialize );
      }

      //! add discrete function to communicator with given unpack operation
      template< class Communicator, class Operation >
      void addToList ( Communicator &comm, const Operation& op)
      {
        rpOp_.addToList( comm, op );
      }

      //! add discrete function to communicator
      template< class Communicator >
      void addToList ( Communicator &comm  )
      {
        rpOp_.addToList( comm );
      }

      //! remove discrete function from communicator
      template< class Communicator >
      void removeFromList ( Communicator &comm )
      {
        rpOp_.removeFromList( comm );
      }

      //! add discrete function to load balancer
      template< class LoadBalancer >
      void addToLoadBalancer ( LoadBalancer& lb )
      {
        rpOp_.addToLoadBalancer( lb );
      }

      void initialize ()
      {
        adaptiveFunction_.assign( discreteFunction_ );
        rpOp_.initialize();
        initialized_ = true ;
      }

      void finalize   ()
      {
        // only finalize if previously initialized
        assert( initialized_ );
        rpOp_.finalize();
        discreteFunction_.assign( adaptiveFunction_ );
        initialized_ = false ;
      }

    protected:
      using BaseType::calcWeight;
      using BaseType::entitiesAreCopies;

    protected:
      DiscreteFunctionType& discreteFunction_;
      AdaptiveDiscreteFunctionType adaptiveFunction_;
      AdaptiveRestrictProlongType  rpOp_;

      bool initialized_;
    };

  } // namespace Fem
} // namespace Dune

#endif // #if HAVE_PETSC

#endif // #ifndef DUNE_FEM_PETSCDISCRETEFUNCTION_HH
