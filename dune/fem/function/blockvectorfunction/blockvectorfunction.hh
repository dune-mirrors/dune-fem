#ifndef DUNE_FEM_BLOCKVECTORFUNCTION_HH
#define DUNE_FEM_BLOCKVECTORFUNCTION_HH

#include <memory>
#include <string>
#include <utility>

#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/function/blockvectorfunction/declaration.hh>
#include <dune/fem/function/blockvectors/defaultblockvectors.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune
{
  namespace Fem
  {
    /** \class DiscreteFunctionTraits
     *  \brief Traits class for a DiscreteFunction
     *
     *  \tparam  DiscreteFunctionSpace   space the discrete function lives in
     *  \tparam  DofVector             implementation class of the block vector
     */
    template< class DiscreteFunctionSpace,
              class Block >
    struct DiscreteFunctionTraits< ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpace, Block > >
      : public DefaultDiscreteFunctionTraits< DiscreteFunctionSpace, ISTLBlockVector< Block > >
    {
      typedef ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpace, Block >  DiscreteFunctionType;
      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };



    template < class DiscreteFunctionSpace, class Block >
    class ISTLBlockVectorDiscreteFunction
    : public DiscreteFunctionDefault< ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpace, Block > >
    {
      typedef ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpace, Block > ThisType;
      typedef DiscreteFunctionDefault< ThisType > BaseType;

    public:
      typedef typename BaseType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
      typedef typename BaseType :: DofVectorType              DofVectorType;
      typedef typename BaseType :: DofType                    DofType;
      typedef typename DofVectorType :: DofContainerType      DofContainerType;
      typedef DofContainerType                                DofStorageType;

      typedef typename BaseType :: ScalarProductType          ScalarProductType;

      using BaseType::assign;

      /** \brief Constructor to use if the vector storing the dofs does not exist yet
       *
       *  \param[in]  name         name of the discrete function
       *  \param[in]  space        space the discrete function lives in
       */
      ISTLBlockVectorDiscreteFunction( const std::string& name,
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
      ISTLBlockVectorDiscreteFunction( const std::string& name,
                                       const DiscreteFunctionSpaceType& space,
                                       const DofContainerType& dofVector )
        : BaseType( name, space ),
          memObject_(),
          dofVector_( const_cast< DofContainerType& > (dofVector) )
      {}

      /** \brief Copy constructor */
      ISTLBlockVectorDiscreteFunction( const ThisType& other )
        : BaseType( "copy of " + other.name(), other.space() ),
          memObject_(),
          dofVector_( allocateDofStorage( other.space() ) )
      {
        assign( other );
      }

      /** \brief Move constructor */
      ISTLBlockVectorDiscreteFunction( ThisType&& other )
        : BaseType( static_cast< BaseType && >( other ) ),
          memObject_( std::move( other.memObject_ ) ),
          dofVector_( std::move( other.dofVector_ ) )
      {}

      ISTLBlockVectorDiscreteFunction () = delete;
      ThisType& operator= ( const ThisType& ) = delete;
      ThisType& operator= ( ThisType&& ) = delete;

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::enableDofCompression() */
      void enableDofCompression ()
      {
        if( memObject_ )
          memObject_->enableDofCompression();
      }

#if HAVE_PETSC
      void assign( const PetscDiscreteFunction< DiscreteFunctionSpaceType >& g )
      {
        g.dofVector().copyTo( dofVector() );
      }
#endif

      //! convenience method for usage with ISTL solvers
      DofContainerType& blockVector() { return dofVector().array(); }

      //! convenience method for usage with ISTL solvers
      const DofContainerType& blockVector() const { return dofVector().array(); }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dofVector() */
      DofVectorType& dofVector() { return dofVector_; }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dofVector() */
      const DofVectorType& dofVector() const { return dofVector_; }

      /** \brief returns ScalarProduct to be used with ISTLInverseOp */
      ScalarProductType& scalarProduct() { return scalarProduct_ ; }

    protected:
      using BaseType :: scalarProduct_;

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
