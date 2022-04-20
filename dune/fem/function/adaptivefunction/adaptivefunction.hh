#ifndef DUNE_FEM_ADAPTIVEFUNCTION_HH
#define DUNE_FEM_ADAPTIVEFUNCTION_HH

#include <memory>
#include <string>
#include <utility>

#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/blockvectors/defaultblockvectors.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/storage/dynamicarray.hh>

namespace Dune
{
  namespace Fem
  {

    template <class DiscreteFunctionSpace>
    class AdaptiveDiscreteFunction;

#if HAVE_PETSC
    template <class DiscreteFunctionSpace>
    class PetscDiscreteFunction;
#endif

    template< typename DiscreteFunctionSpace >
    struct DiscreteFunctionTraits< AdaptiveDiscreteFunction< DiscreteFunctionSpace > >
      : public DefaultDiscreteFunctionTraits< DiscreteFunctionSpace,
         SimpleBlockVector< StaticArray< typename DiscreteFunctionSpace::RangeFieldType > , DiscreteFunctionSpace::localBlockSize > >
    {
      typedef AdaptiveDiscreteFunction< DiscreteFunctionSpace > DiscreteFunctionType;
      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };



    //! @ingroup AdaptiveDFunction
    //! An adaptive discrete function
    //! This class is comparable to DFAdapt, except that it provides a
    //! specialisation for CombinedSpace objects which provides enriched
    //! functionality (access to subfunctions) and runtime optimisations
    template <class DiscreteFunctionSpace>
    class AdaptiveDiscreteFunction
    : public DiscreteFunctionDefault< AdaptiveDiscreteFunction< DiscreteFunctionSpace > >
    {
      typedef AdaptiveDiscreteFunction< DiscreteFunctionSpace > ThisType;
      typedef DiscreteFunctionDefault< ThisType > BaseType;

    public:
      typedef typename BaseType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
      typedef typename BaseType :: DofVectorType              DofVectorType;
      typedef typename BaseType :: DofType                    DofType;

      // generic assign method
      using BaseType::assign;
      using BaseType::name;

      typedef MutableBlockVector< DynamicArray< DofType >, DiscreteFunctionSpaceType::localBlockSize > MutableDofVectorType;

      /** \brief Constructor to use if the vector storing the dofs does not exist yet
       *
       *  \param[in]  name         name of the discrete function
       *  \param[in]  space        space the discrete function lives in
       */
      AdaptiveDiscreteFunction( const std::string& name,
                                const DiscreteFunctionSpaceType& space )
        : BaseType( name, space ),
          memObject_(),
          dofVector_( allocateDofStorage( space ) )
      {}

      /** \brief Constructor to use if the vector storing the dofs already exists
       *
       *  \param[in]  name         name of the discrete function
       *  \param[in]  space        space the discrete function lives in
       *  \param[in]  data         pointer to data
       */
      AdaptiveDiscreteFunction( const std::string& name,
                                const DiscreteFunctionSpaceType& space,
                                const DofType* data )
        : BaseType( name, space ),
          memObject_(),
          dofVector_( allocateDofStorageWrapper( space.blockMapper().size() * DofVectorType::blockSize, data ) )
      {}

      /** \brief Constructor to use if the vector storing the dofs already exists
       *
       *  \param[in]  name         name of the discrete function
       *  \param[in]  space        space the discrete function lives in
       *  \param[in]  dofVector    reference to the dof vector
       */
      AdaptiveDiscreteFunction( const std::string& name,
                                const DiscreteFunctionSpaceType& space,
                                DofVectorType& dofVector )
        : BaseType( name, space ),
          memObject_(),
          dofVector_( dofVector )
      {}

      /** \brief Copy constructor */
      AdaptiveDiscreteFunction( const ThisType& other )
        : BaseType( "copy of " + other.name(), other.space() ),
          memObject_(),
          dofVector_( allocateDofStorage( other.space() ) )
      {
        assign( other );
      }

      /** \brief Move constructor */
      AdaptiveDiscreteFunction( ThisType&& other )
        : BaseType( static_cast< BaseType && >( other ) ),
          memObject_( std::move( other.memObject_ ) ),
          dofVector_( other.dofVector_ )
      {}

      AdaptiveDiscreteFunction () = delete;
      ThisType& operator= ( const ThisType& ) = delete;
      ThisType& operator= ( ThisType&& ) = delete;

#if HAVE_PETSC
      void assign( const PetscDiscreteFunction< DiscreteFunctionSpaceType >& g )
      {
        g.dofVector().copyTo( dofVector() );
      }
#endif

      DofType* leakPointer() { return dofVector().data(); }
      const DofType* leakPointer() const { return dofVector().data(); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dofVector() */
      DofVectorType& dofVector() { return dofVector_; }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dofVector() */
      const DofVectorType& dofVector() const { return dofVector_; }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::enableDofCompression() */
      void enableDofCompression ()
      {
        if( memObject_ )
          memObject_->enableDofCompression();
      }

    protected:

      //! wrapper class to create fake DofStorage from DofType*
      class DofStorageWrapper : public DofStorageInterface
      {
        StaticArray< DofType > array_;
        DofVectorType dofVector_;
        typedef typename DofStorageInterface::SizeType SizeType;

      public:
        DofStorageWrapper ( const SizeType size,
                            const DofType *v )
          : array_( size, const_cast< DofType* >(v) ),
            dofVector_( array_ )
        {}

        //! return array
        DofVectorType &getArray () { return dofVector_; }

        //! do nothing here since we are using StaticArray
        void enableDofCompression () {}

        //! return array's size
        SizeType size () const { return dofVector_.size(); }
      };

    protected:
      // allocate unmanaged dof storage
      DofVectorType&
      allocateDofStorageWrapper ( const size_t size,
                                  const DofType *v )
      {
        DofStorageWrapper *dsw = new DofStorageWrapper( size, v );
        assert( dsw );

        // save pointer to object
        memObject_.reset( dsw );
        // return array
        return dsw->getArray();
      }


      // allocate managed dof storage
      DofVectorType& allocateDofStorage ( const DiscreteFunctionSpaceType &space )
      {
        // create memory object
        std::pair< DofStorageInterface*, DofVectorType* > memPair
          = allocateManagedDofStorage( space.gridPart().grid(), space.blockMapper(),
                                       (MutableDofVectorType *) 0 );

        // save pointer
        memObject_.reset( memPair.first );
        return *(memPair.second);
      }

      std::unique_ptr< DofStorageInterface > memObject_;
      DofVectorType& dofVector_;
    };

  } // end namespace Fem
} // end namespace Dune

#endif // #ifndef DUNE_FEM_ADAPTIVEFUNCTION_HH
