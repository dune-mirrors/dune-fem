#ifndef DUNE_FEM_BLOCKVECTORFUNCTION_HH
#define DUNE_FEM_BLOCKVECTORFUNCTION_HH

//- system includes
#include <fstream>
#include <iostream>

//- Dune inlcudes
#include <dune/common/exceptions.hh>
#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/space/common/dofmanager.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bvector.hh>
#else
#include <dune/fem/space/common/arrays.hh>
#endif

#include <dune/fem/function/blockvectors/defaultblockvectors.hh>

#include <dune/fem/common/referencevector.hh>
#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/localfunction/mutable.hh>

namespace Dune
{

  namespace Fem
  {
    //**********************************************************************
    //! @ingroup BlockVectorDFunction
    //  --ISTLBlockVectorDiscreteFunction
    //
    //! this is one special implementation of a discrete function using an
    //! array for storing the dofs.
    //!
    //**********************************************************************
    template <class DiscreteFunctionSpace>
    class ISTLBlockVectorDiscreteFunction;

    /** \class DiscreteFunctionTraits
     *  \brief Traits class for a DiscreteFunction
     *
     *  \tparam  DiscreteFunctionSpace   space the discrete function lives in
     *  \tparam  DofVector             implementation class of the block vector
     */
    template< typename DiscreteFunctionSpace >
    struct DiscreteFunctionTraits< ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpace > >
      : public DefaultDiscreteFunctionTraits< DiscreteFunctionSpace,
            ISTLBlockVector< Dune::FieldVector< typename DiscreteFunctionSpace::RangeFieldType, DiscreteFunctionSpace::localBlockSize > > >
    {
      typedef ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpace >  DiscreteFunctionType;
      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };


    template <class DiscreteFunctionSpace>
    class ISTLBlockVectorDiscreteFunction
    : public DiscreteFunctionDefault< ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpace > >
    {
      typedef ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpace > ThisType;
      typedef DiscreteFunctionDefault< ThisType > BaseType;

    public:
      typedef typename BaseType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
      typedef typename BaseType :: DofVectorType              DofVectorType;
      typedef typename BaseType :: DofType                    DofType;
      typedef typename DofVectorType :: DofContainerType      DofContainerType;
      typedef DofContainerType                                DofStorageType;

      using BaseType::assign;

      ISTLBlockVectorDiscreteFunction( const std::string &name,
                                       const DiscreteFunctionSpaceType &space )
        : BaseType( name, space ),
          memObject_(),
          dofVector_( allocateDofStorage( space ) )
      {
      }

      ISTLBlockVectorDiscreteFunction( const std::string &name,
                                       const DiscreteFunctionSpaceType &space,
                                       const DofContainerType& dofContainer )
        : BaseType( name, space ),
          memObject_(),
          dofVector_( const_cast< DofContainerType& > (dofContainer) )
      {
      }

      ISTLBlockVectorDiscreteFunction( const ISTLBlockVectorDiscreteFunction& other )
        : BaseType( "copy of " + other.name(), other.space() ),
          memObject_(),
          dofVector_( allocateDofStorage( other.space() ) )
      {
        assign( other );
      }

      void enableDofCompression ()
      {
        if( memObject_ )
          memObject_->enableDofCompression();
      }

      //! convenience method for usage with ISTL solvers
      DofContainerType& blockVector() { return dofVector().array(); }
      //! convenience method for usage with ISTL solvers
      const DofContainerType& blockVector() const { return dofVector().array(); }

      DofVectorType& dofVector() { return dofVector_; }
      const DofVectorType& dofVector() const { return dofVector_; }

    protected:
      // allocate managed dof storage
      DofContainerType& allocateDofStorage ( const DiscreteFunctionSpaceType &space )
      {
        std::string name("deprecated");
        // create memory object
        std::pair< DofStorageInterface*, DofContainerType* > memPair
          = allocateManagedDofStorage( space.gridPart().grid(), space.blockMapper(),
                                       name, (DofContainerType *) 0 );

        // save pointer
        memObject_.reset( memPair.first );
        return *(memPair.second);
      }

      // pointer to allocated DofVector
      std::unique_ptr< DofStorageInterface > memObject_;

      // DofVector object holds pointer to dof container
      DofVectorType dofVector_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BLOCKVECTORFUNCTION_HH
