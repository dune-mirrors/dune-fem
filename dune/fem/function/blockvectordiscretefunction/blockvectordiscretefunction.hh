// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_BLOCKVECTORDISCRETEFUNCTION_HH
#define DUNE_FEM_BLOCKVECTORDISCRETEFUNCTION_HH

#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/common/referencevector.hh>
#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/function/localfunction/mutable.hh>

#include <dune/fem/misc/threads/threadmanager.hh>
#include <dune/fem/storage/envelope.hh>

namespace Dune
{

  namespace Fem
  {

    // forward declaration
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

      /*
         I didn't implement these methods of DiscreteFunctionDefault (deliberately):
          void print ( std :: ostream &out ) const;
       */
      typedef BlockVectorDiscreteFunction< DiscreteFunctionSpace, BlockVector >   ThisType;
      typedef DiscreteFunctionDefault< ThisType > BaseType;

      typedef ParallelScalarProduct< ThisType >                                   ScalarProductType;

    public:
      // ==================== Types

      //! type for the discrete function space this function lives in
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
      //! type for the class which implements the block vector
      typedef BlockVector BlockVectorType;
      //! type for the class which implements the block vector (which is the dof vector)
      typedef BlockVectorType DofVectorType;

      // methods from DiscreteFunctionDefault
      using BaseType::assign;

      /** \brief Constructor to use if the vector storing the dofs (which is a block vector) already exists
       *
       *  \param[in]  name         name of the discrete function
       *  \param[in]  dfSpace      space the discrete function lives in
       *  \param[in]  blockVector  reference to the blockVector
       */
      BlockVectorDiscreteFunction ( const std::string &name,
                                    const DiscreteFunctionSpaceType &dfSpace,
                                    DofVectorType &dofVector )
        : BaseType( name, dfSpace ),
          memObject_(),
          dofVector_( dofVector )
      {}

      /** \brief Constructor to use if the vector storing the dofs does not exist yet
       *
       *  \param[in]  name         name of the discrete function
       *  \param[in]  dfSpace      space the discrete function lives in
       */
      BlockVectorDiscreteFunction ( const std::string &name,
                                    const DiscreteFunctionSpaceType &dfSpace )
        : BaseType( name, dfSpace ),
          memObject_(),
          dofVector_( allocateDofStorage( dfSpace ) )
      {}


      /** \brief Copy constructor
       */
      BlockVectorDiscreteFunction ( const ThisType &other )
        : BaseType( "copy of "+other.name(), other.space() ),
          memObject_(),
          dofVector_( allocateDofStorage( other.space() ) )
      {
        // copy dof vector content
        assign( other );
      }

    private:
      // an empty constructor would not make sense for a discrete function
      BlockVectorDiscreteFunction ();

      // TODO: un-disallow this??
      ThisType &operator= ( const ThisType &other );

    public:
      /** \brief Obtain constant reference to the dof vector
       *
       *  \return Constant reference to the block vector
       */
      const DofVectorType &dofVector () const
      {
        return dofVector_;
      }

      /** \brief Obtain reference to the dof vector
       *
       *  \return Reference to the block vector
       */
      DofVectorType &dofVector ()
      {
        return dofVector_;
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::enableDofCompression()
       */
      void enableDofCompression ()
      {
        if( memObject_ )
          memObject_->enableDofCompression();
      }

    protected:
      DofVectorType& allocateDofStorage( const DiscreteFunctionSpaceType& space )
      {
        std::string name("deprecated");
        std::pair< DofStorageInterface*, DofVectorType* > memPair(
          allocateManagedDofStorage< DofVectorType >( space.gridPart().grid(), space.blockMapper(), name ) );

        memObject_.reset( memPair.first );
        return *memPair.second;
      }

      /*
       * ============================== data fields ====================
       */
      std::unique_ptr< DofStorageInterface > memObject_;
      DofVectorType& dofVector_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BLOCKVECTORDISCRETEFUNCTION_HH
