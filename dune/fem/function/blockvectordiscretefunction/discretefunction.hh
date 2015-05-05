// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_NORMALDISCRETEFUNCTION_HH
#define DUNE_FEM_NORMALDISCRETEFUNCTION_HH

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
    template< typename DiscreteFunctionSpace, typename DofVector >
    class DiscreteFunction;

    /** \class DiscreteFunctionTraits
     *  \brief Traits class for a DiscreteFunction
     *
     *  \tparam  DiscreteFunctionSpace   space the discrete function lives in
     *  \tparam  DofVector             implementation class of the block vector
     */
    template< typename DiscreteFunctionSpace, typename DofVector >
    struct DiscreteFunctionTraits< DiscreteFunction< DiscreteFunctionSpace, DofVector > >
    {
      typedef DofVector DofVectorType;

      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
      typedef DiscreteFunction< DiscreteFunctionSpaceType, DofVectorType >   DiscreteFunctionType;

      typedef typename DofVectorType::IteratorType      DofIteratorType;
      typedef typename DofVectorType::ConstIteratorType ConstDofIteratorType;
      typedef typename DofVectorType::DofBlockType      DofBlockType;
      typedef typename DofVectorType::ConstDofBlockType ConstDofBlockType;
      typedef typename DofVectorType::DofContainerType  DofContainerType;

      typedef Fem::Envelope< DofBlockType >             DofBlockPtrType;
      typedef Fem::Envelope< ConstDofBlockType >        ConstDofBlockPtrType;

      typedef typename DiscreteFunctionSpaceType::BlockMapperType MapperType;
      typedef typename DofVectorType::FieldType DofType;

      typedef ThreadSafeValue< UninitializedObjectStack > LocalDofVectorStackType;
      typedef StackAllocator< DofType, LocalDofVectorStackType* > LocalDofVectorAllocatorType;
      typedef DynamicReferenceVector< DofType, LocalDofVectorAllocatorType > LocalDofVectorType;

      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };

    /** \class DiscreteFunctionTraits
     *  \brief A discrete function which uses block vectors for its dof storage and dof calculations
     *
     *  \tparam  DiscreteFunctionSpace   space the discrete function lives in
     *  \tparam  DofVector             implementation class of the block vector
     */
    template< typename DiscreteFunctionSpace, typename DofVector >
    class DiscreteFunction
      : public DiscreteFunctionDefault< DiscreteFunction< DiscreteFunctionSpace, DofVector > >
    {

      /*
         I didn't implement these methods of DiscreteFunctionDefault (deliberately):
          void print ( std :: ostream &out ) const;
       */
      typedef DiscreteFunction< DiscreteFunctionSpace, DofVector >   ThisType;
      typedef DiscreteFunctionDefault< DiscreteFunction< DiscreteFunctionSpace, DofVector > > BaseType;

      typedef ParallelScalarProduct< ThisType >                                   ScalarProductType;

    public:
      // ==================== Types

      //! the traits of ThisType
      typedef DiscreteFunctionTraits< ThisType > Traits;
      //! type for the discrete function space this function lives in
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
      //! type for the class which implements the block vector
      typedef DofVector DofVectorType;
      //! type for the mapper which maps local to global indices
      typedef typename Traits::MapperType MapperType;
      //! type for the mapper which maps local to global indices
      typedef typename Traits::DofContainerType DofContainerType;

      //! type of the fields the dofs live in
      typedef typename BaseType::DofType DofType;
      //! pointer to a block of dofs
      typedef typename BaseType::DofBlockPtrType DofBlockPtrType;
      //! pointer to a block of dofs, const version
      typedef typename BaseType::ConstDofBlockPtrType ConstDofBlockPtrType;
      //! type of a block of dofs
      typedef typename BaseType::DofBlockType DofBlockType;
      //! iterator type for iterate over the dofs
      typedef typename BaseType::DofIteratorType DofIteratorType;
      //! iterator type for iterate over the dofs, const verision
      typedef typename BaseType::ConstDofIteratorType ConstDofIteratorType;
      //! type of the range field
      typedef typename BaseType::RangeFieldType RangeFieldType;
      //! type of the domain field
      typedef typename BaseType::DomainFieldType DomainFieldType;
      //! type of the local functions
      typedef typename BaseType::LocalFunctionType LocalFunctionType;
      //! type of LocalDofVector
      typedef typename BaseType::LocalDofVectorType LocalDofVectorType;
      //! type of LocalDofVector StackAllocator
      typedef typename BaseType::LocalDofVectorAllocatorType LocalDofVectorAllocatorType;
      //! type of the discrete functions's domain
      typedef typename BaseType::DomainType DomainType;
      //! type of the discrete functions's range
      typedef typename BaseType::RangeType RangeType;
      //! size type of the block vector
      typedef typename DofVectorType::SizeType SizeType;

      // methods from DiscreteFunctionDefault
      using BaseType::space;
      using BaseType::name;
      using BaseType::assign;
      using BaseType::operator+=;
      using BaseType::operator-=;

      //! size of the dof blocks
      enum { blockSize = DiscreteFunctionSpaceType::localBlockSize };

      /** \brief Constructor to use if the vector storing the dofs (which is a block vector) already exists
       *
       *  \param[in]  name         name of the discrete function
       *  \param[in]  dfSpace      space the discrete function lives in
       *  \param[in]  blockVector  reference to the blockVector
       */
      DiscreteFunction ( const std::string &name,
                         const DiscreteFunctionSpaceType &dfSpace,
                         DofVectorType &dofVector )
        : BaseType( name, dfSpace, LocalDofVectorAllocatorType( &ldvStack_ ) ),
          ldvStack_( std::max( std::max( sizeof( DofType ), sizeof( DofType* ) ),
                               sizeof(typename LocalDofVectorType::value_type) ) // for PetscDiscreteFunction
                        * space().blockMapper().maxNumDofs() * DiscreteFunctionSpaceType::localBlockSize ),
          dofVector_( dofVector )
      {}

    private:

      // an empty constructor would not make sense for a discrete function
      DiscreteFunction ();

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

      /** \brief Obtain constant reference to the dof vector
       *
       *  \return Constant reference to the block vector
       */
      const DofContainerType &dofContainer () const
      {
        return dofVector().array();
      }

      /** \brief Obtain reference to the dof vector
       *
       *  \return Reference to the block vector
       */
      DofContainerType &dofContainer ()
      {
        return dofVector().array();
      }

      /** \brief Return the number of blocks in the block vector
       *
       *  \return Number of block in the block vector
       */
      SizeType size () const { return dofVector().size() * blockSize; }

    protected:
      /*
       * ============================== data fields ====================
       */
      typename Traits :: LocalDofVectorStackType ldvStack_;
      DofVectorType& dofVector_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BLOCKVECTORDISCRETEFUNCTION_HH
