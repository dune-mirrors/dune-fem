// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_PETSCDISCRETEFUNCTION_HH
#define DUNE_FEM_PETSCDISCRETEFUNCTION_HH

#include <string>
#include <algorithm>
#include <map>
#include <vector>
#include <iterator>
#include <utility>
#include <memory>

#if HAVE_PETSC

#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/misc/petsc/petscdofmappings.hh>
#include <dune/fem/misc/petsc/petscvector.hh>
#include <dune/fem/misc/petsc/petscslavedofprovider.hh>

#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>

#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/functor.hh>
#include <dune/fem/function/localfunction/mutable.hh>

namespace Dune
{

  namespace Fem
  {


    /* ==============================
     * forward declarations
     */
    template< class DiscreteFunctionSpace > class PetscDiscreteFunction;

    // AssignVectorReference
    // ---------------------

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

    // Internal Forward Declaration
    //-----------------------------

    template <class DiscreteFunctionSpace>
    class PetscDiscreteFunction;

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
    //----------------------

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

      PetscDiscreteFunction( const std::string &name,
                             const DiscreteFunctionSpaceType &space )
        : BaseType( name, space ),
          memObject_(),
          dofVector_( allocateDofStorage( space ) )
      {
      }

      PetscDiscreteFunction( const std::string &name,
                             const DiscreteFunctionSpaceType &space,
                             DofVectorType& dofVector )
        : BaseType( name, space ),
          memObject_(),
          dofVector_( dofVector )
      {
      }

      PetscDiscreteFunction( const PetscDiscreteFunction& other )
        : BaseType( "copy of " + other.name(), other.space() ),
          memObject_(),
          dofVector_( allocateDofStorage( other.space() ) )
      {
        assign( other );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::enableDofCompression()
       */
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
