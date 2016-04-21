#ifndef DUNE_FEM_ADAPTIVEFUNCTION_HH
#define DUNE_FEM_ADAPTIVEFUNCTION_HH

//- System includes
#include <memory>
#include <string>
#include <utility>

//- Dune includes
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/blockvectors/defaultblockvectors.hh>

//- Local includes
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/common/referencevector.hh>
#include <dune/fem/common/stackallocator.hh>

namespace Dune
{

  namespace Fem
  {

    //! @ingroup AdaptiveDFunction
    //! An adaptive discrete function
    //! This class is comparable to DFAdapt, except that it provides a
    //! specialisation for CombinedSpace objects which provides enriched
    //! functionality (access to subfunctions) and runtime optimisations
    template <class DiscreteFunctionSpace>
    class AdaptiveDiscreteFunction;

    template< typename DiscreteFunctionSpace >
    struct DiscreteFunctionTraits< AdaptiveDiscreteFunction< DiscreteFunctionSpace > >
      : public DefaultDiscreteFunctionTraits< DiscreteFunctionSpace,
         SimpleBlockVector< StaticArray< typename DiscreteFunctionSpace::RangeFieldType > , DiscreteFunctionSpace::localBlockSize > >
    {
      typedef AdaptiveDiscreteFunction< DiscreteFunctionSpace > DiscreteFunctionType;
      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };

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

      using BaseType::assign;

      typedef MutableBlockVector< MutableArray< DofType >, DiscreteFunctionSpaceType::localBlockSize > MutableDofVectorType;

      AdaptiveDiscreteFunction( const std::string &name,
                                const DiscreteFunctionSpaceType &space )
        : BaseType( name, space ),
          memObject_(),
          dofVector_( allocateDofStorage( space ) )
      {
      }

      AdaptiveDiscreteFunction( const std::string &name,
                                const DiscreteFunctionSpaceType &space,
                                const DofType* data )
        : BaseType( name, space ),
          memObject_(),
          dofVector_( allocateDofStorageWrapper( space.blockMapper().size() * DofVectorType::blockSize, data ) )
      {
      }

      AdaptiveDiscreteFunction( const std::string &name,
                                const DiscreteFunctionSpaceType &space,
                                DofVectorType& dofVector )
        : BaseType( name, space ),
          memObject_(),
          dofVector_( dofVector )
      {
      }

      AdaptiveDiscreteFunction( const AdaptiveDiscreteFunction& other )
        : BaseType( "copy of " + other.name(), other.space() ),
          memObject_(),
          dofVector_( allocateDofStorage( other.space() ) )
      {
        assign( other );
      }

      DofType* leakPointer() { return dofVector().data(); }
      const DofType* leakPointer() const { return dofVector().data(); }

      DofVectorType& dofVector() { return dofVector_; }
      const DofVectorType& dofVector() const { return dofVector_; }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::enableDofCompression()
       */
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
        typedef typename DofVectorType :: SizeType SizeType;

        std::string name_;

      public:
        DofStorageWrapper ( const SizeType size,
                            const DofType *v )
          : array_( size, const_cast< DofType* >(v) ),
            dofVector_( array_ ),
            name_("deprecated")
        {}

        const std::string& name () const { return name_; }

        //! return array
        DofVectorType &getArray () { return dofVector_; }

        //! do nothing here since we are using StaticArray
        void enableDofCompression () {}

        //! return array's size
        int size () const { return dofVector_.size(); }
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
        std::string name("deprecated");
        // create memory object
        std::pair< DofStorageInterface*, DofVectorType* > memPair
          = allocateManagedDofStorage( space.gridPart().grid(), space.blockMapper(),
                                       name, (MutableDofVectorType *) 0 );

        // save pointer
        memObject_.reset( memPair.first );
        return *(memPair.second);
      }

      // pointer to allocated DofVector
      std::unique_ptr< DofStorageInterface > memObject_;

      DofVectorType& dofVector_;
    };

  } // end namespace Fem

} // end namespace Dune
#endif // #ifndef DUNE_FEM_ADAPTIVEFUNCTION_HH
