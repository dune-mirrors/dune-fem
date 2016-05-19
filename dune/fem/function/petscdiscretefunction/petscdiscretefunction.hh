#ifndef DUNE_FEM_PETSCDISCRETEFUNCTION_HH
#define DUNE_FEM_PETSCDISCRETEFUNCTION_HH

#include <memory>
#include <string>
#include <utility>

#if HAVE_PETSC

#include <dune/common/dynvector.hh>
#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/functor.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/misc/petsc/petscdofmappings.hh>
#include <dune/fem/misc/petsc/petscvector.hh>
#include <dune/fem/misc/petsc/petscslavedofprovider.hh>

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

      void operator() ( const std::size_t index, DofProxy value )
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



    template <class DiscreteFunctionSpace>
    class PetscDiscreteFunction
    : public DiscreteFunctionDefault< PetscDiscreteFunction< DiscreteFunctionSpace > >
    {
      typedef PetscDiscreteFunction< DiscreteFunctionSpace > ThisType;
      typedef DiscreteFunctionDefault< ThisType > BaseType;

    public:
      typedef typename BaseType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
      typedef typename BaseType :: DofVectorType              DofVectorType;

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

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::communicate() */
      void communicate ()
      {
        dofVector().communicateNow();
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dofVector() */
      DofVectorType& dofVector() { return dofVector_; }
      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dofVector() */
      const DofVectorType& dofVector() const { return dofVector_; }

      /** \brief obtain a constand pointer to the underlying PETSc Vec */
      const Vec* petscVec () const { return dofVector().getVector(); }

      /** \brief obtain a pointer to the underlying PETSc Vec */
      Vec* petscVec () { return dofVector().getVector(); }

    protected:
      typedef typename DiscreteFunctionSpaceType :: BlockMapperType  BlockMapperType;
      typedef PetscManagedDofStorage< DiscreteFunctionSpaceType, BlockMapperType > PetscManagedDofStorageType;

      // allocate managed dof storage
      DofVectorType& allocateDofStorage ( const DiscreteFunctionSpaceType &space )
      {
        std::string name("deprecated");

        memObject_.reset( new PetscManagedDofStorageType( space, space.blockMapper(), name ) );

        return memObject_->getArray();
      }

      // pointer to allocated DofVector
      std::unique_ptr< PetscManagedDofStorageType > memObject_;

      // dof vector impl
      DofVectorType& dofVector_;
    };

  } // namespace Fem
} // namespace Dune

#endif // #if HAVE_PETSC

#endif // #ifndef DUNE_FEM_PETSCDISCRETEFUNCTION_HH
