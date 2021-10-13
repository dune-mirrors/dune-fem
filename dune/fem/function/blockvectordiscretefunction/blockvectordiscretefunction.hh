#ifndef DUNE_FEM_BLOCKVECTORDISCRETEFUNCTION_HH
#define DUNE_FEM_BLOCKVECTORDISCRETEFUNCTION_HH

#include <memory>
#include <string>
#include <utility>

#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/function/blockvectors/referenceblockvector.hh>
#include <dune/fem/storage/envelope.hh>

namespace Dune
{
  namespace Fem
  {

    template< typename DiscreteFunctionSpace, typename BlockVector >
    class BlockVectorDiscreteFunction;



    /** \class IsBlockVectorDiscreteFunction
     *  \brief Tag for discrete functions using block vectors
     *
     *  A discrete function using block vectors for its dof storage and calculations should inherit from
     *  this struct. For example, Dune::Fem::MatrixOperator recognizes discrete functions with block vectors
     *  only as such if they inherit from this tag. If they do, their method .dofVector() is used (which is
     *  a block vector). If they do not, this indicates that they don't use block vectors and thus provide
     *  a .leakPointer() method -  which is used by Dune::Fem::MatrixOperator in this case.
     *
     */
    struct IsBlockVectorDiscreteFunction {};



    /** \class BlockVectorDiscreteFunctionTraits
     *  \brief Traits class for a BlockVectorDiscreteFunction
     *
     *  \tparam  DiscreteFunctionSpace   space the discrete function lives in
     *  \tparam  BlockVector             implementation class of the block vector
     */
    template< typename DiscreteFunctionSpace, typename BlockVector >
    struct DiscreteFunctionTraits< BlockVectorDiscreteFunction< DiscreteFunctionSpace, BlockVector > >
      : public DefaultDiscreteFunctionTraits< DiscreteFunctionSpace, BlockVector >
    {
      typedef BlockVectorDiscreteFunction< DiscreteFunctionSpace, BlockVector >   DiscreteFunctionType;
      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };




    /** \class BlockVectorDiscreteFunctionTraits
     *  \brief A discrete function which uses block vectors for its dof storage and dof calculations
     *
     *  \tparam  DiscreteFunctionSpace   space the discrete function lives in
     *  \tparam  BlockVector             implementation class of the block vector
     */
    template< typename DiscreteFunctionSpace, typename BlockVector >
    class BlockVectorDiscreteFunction
      : public DiscreteFunctionDefault< BlockVectorDiscreteFunction< DiscreteFunctionSpace, BlockVector > >,
        public IsBlockVectorDiscreteFunction
    {
      typedef BlockVectorDiscreteFunction< DiscreteFunctionSpace, BlockVector > ThisType;
      typedef DiscreteFunctionDefault< ThisType > BaseType;
      typedef ParallelScalarProduct< ThisType > ScalarProductType;

    public:
      //! type for the discrete function space this function lives in
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
      //! type for the class which implements the block vector
      typedef BlockVector BlockVectorType;
      //! type for the class which implements the block vector (which is the dof vector)
      typedef BlockVectorType DofVectorType;

      using BaseType::assign;
      using BaseType::name;

      /** \brief Constructor to use if the vector storing the dofs (which is a block vector) already exists
       *
       *  \param[in]  name         name of the discrete function
       *  \param[in]  space        space the discrete function lives in
       *  \param[in]  blockVector  reference to the blockVector
       */
      BlockVectorDiscreteFunction ( const std::string& name,
                                    const DiscreteFunctionSpaceType& space,
                                    DofVectorType& dofVector )
        : BaseType( name, space ),
          memObject_(),
          dofVector_( dofVector )
      {}

      /** \brief Constructor to use if the vector storing the dofs does not exist yet
       *
       *  \param[in]  name         name of the discrete function
       *  \param[in]  space        space the discrete function lives in
       */
      BlockVectorDiscreteFunction ( const std::string& name,
                                    const DiscreteFunctionSpaceType& space )
        : BaseType( name, space ),
          memObject_(),
          dofVector_( allocateDofStorage( space ) )
      {}

      /** \brief Copy constructor */
      BlockVectorDiscreteFunction ( const ThisType& other )
        : BaseType( "copy of "+other.name(), other.space() ),
          memObject_(),
          dofVector_( allocateDofStorage( other.space() ) )
      {
        assign( other );
      }

      /** \brief Move constructor */
      BlockVectorDiscreteFunction( ThisType&& other )
        : BaseType( static_cast< BaseType && >( other ) ),
          memObject_( std::move( other.memObject_ ) ),
          dofVector_( other.dofVector_ )
      {}

      BlockVectorDiscreteFunction () = delete;
      ThisType& operator= ( const ThisType& ) = delete;
      ThisType& operator= ( ThisType&& ) = delete;

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dofVector() */
      const DofVectorType &dofVector () const
      {
        return dofVector_;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dofVector() */
      DofVectorType &dofVector ()
      {
        return dofVector_;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::enableDofCompression() */
      void enableDofCompression ()
      {
        if( memObject_ )
          memObject_->enableDofCompression();
      }

    protected:
      DofVectorType& allocateDofStorage( const DiscreteFunctionSpaceType& space )
      {
        std::pair< DofStorageInterface*, DofVectorType* > memPair(
          allocateManagedDofStorage< DofVectorType >( space.gridPart().grid(), space.blockMapper() ) );

        memObject_.reset( memPair.first );
        return *memPair.second;
      }

      std::unique_ptr< DofStorageInterface > memObject_;
      DofVectorType& dofVector_;
    };

  } // namespace Fem
} // namespace Dune

#endif // #ifndef DUNE_FEM_BLOCKVECTORDISCRETEFUNCTION_HH
