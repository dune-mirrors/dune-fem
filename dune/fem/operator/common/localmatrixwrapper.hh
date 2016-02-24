#ifndef DUNE_FEM_LOCALMATRIXWRAPPER_HH
#define DUNE_FEM_LOCALMATRIXWRAPPER_HH

#include <dune/fem/operator/common/localmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class LocalMatrixStack >
    class LocalMatrixWrapper;



    // LocalMatrixWrapperTraits
    // ------------------------

    template< class LocalMatrixStack >
    struct LocalMatrixWrapperTraits
    {
      typedef LocalMatrixStack LocalMatrixStackType;

      typedef typename LocalMatrixStack::ObjectType WrappedLocalMatrixType;

      typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;

      typedef typename WrappedLocalMatrixType::RangeFieldType RangeFieldType;

      typedef typename WrappedLocalMatrixType::DomainSpaceType DomainSpaceType;
      typedef typename WrappedLocalMatrixType::RangeSpaceType RangeSpaceType;

      typedef typename WrappedLocalMatrixType::LittleBlockType LittleBlockType;
    };



    // LocalMatrixWrapper
    // ------------------

    template< class LocalMatrixStack >
    class LocalMatrixWrapper
    : public LocalMatrixInterface< LocalMatrixWrapperTraits< LocalMatrixStack > >
    {
      typedef LocalMatrixWrapper< LocalMatrixStack > ThisType;
      typedef LocalMatrixInterface< LocalMatrixWrapperTraits< LocalMatrixStack > > BaseType;

    public:
      //! type of the local matrix stack
      typedef LocalMatrixStack LocalMatrixStackType;

      //! type of the traits
      typedef LocalMatrixWrapperTraits< LocalMatrixStackType > Traits;

      //! type of the wrapped local matrix
      typedef typename Traits::WrappedLocalMatrixType WrappedLocalMatrixType;

      typedef typename Traits::RangeFieldType RangeFieldType;

      typedef typename BaseType::DomainSpaceType DomainSpaceType;
      typedef typename BaseType::RangeSpaceType RangeSpaceType;

      typedef typename BaseType::DomainBasisFunctionSetType DomainBasisFunctionSetType;
      typedef typename BaseType::RangeBasisFunctionSetType RangeBasisFunctionSetType;

      typedef typename BaseType::DomainEntityType DomainEntityType;
      typedef typename BaseType::RangeEntityType RangeEntityType;

    private:
      typedef typename LocalMatrixStackType::PointerType WrappedLocalMatrixPtrType;

      // ObjectPointer to the actual local matrix
      // (the pointer is required to keep the reference alive)
      WrappedLocalMatrixPtrType localMatrixPtr_;

      // reference to the actual local matrix
      WrappedLocalMatrixType &localMatrix_;

    public:
      //! constructor creating an uninitialized local matrix
      explicit LocalMatrixWrapper ( LocalMatrixStackType &stack )
      : localMatrixPtr_( stack.getObject() ),
        localMatrix_( *localMatrixPtr_ )
      {}

      //! constructor initializing the wrapped local matrix
      template< class DomainEntityType, class RangeEntityType >
      LocalMatrixWrapper( LocalMatrixStackType &stack,
                                 const DomainEntityType &domainEntity,
                                 const RangeEntityType &rangeEntity )
      : localMatrixPtr_( stack.getObject() ),
        localMatrix_( *localMatrixPtr_ )
      {
        // initialize the wrapped local matrix with the entities
        localMatrix().init( domainEntity, rangeEntity );
      }

      /** \brief copy constructor
       *
       *  \param[in]  other  LocalMatrixWrapper to copy
       */
      LocalMatrixWrapper ( const ThisType &other )
      : localMatrixPtr_( other.localMatrixPtr_ ),
        localMatrix_( *localMatrixPtr_ )
      {}

      /** \brief destructor */
      ~LocalMatrixWrapper ( )
      {
        // call finalize on local matrix implementation
        // (e.g. needed for PETSc to add values to the real matrix)
        localMatrix().finalize();
      }

      ThisType& operator= ( const ThisType& ) = delete;

      /** \copydoc Dune::Fem::LocalMatrixInterface::init */
      void init ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity )
      {
        localMatrix().init( domainEntity, rangeEntity );
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::add */
      void add ( int localRow, int localCol, const RangeFieldType &value )
      {
        localMatrix().add( localRow, localCol, value );
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::set */
      void set ( int localRow, int localCol, const RangeFieldType &value )
      {
        localMatrix().set( localRow, localCol, value );
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::clearRow */
      void clearRow ( const int localRow )
      {
        localMatrix().clearRow( localRow );
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::clearRow */
      void clearCol ( const int localCol )
      {
        localMatrix().clearCol( localCol );
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::get */
      const RangeFieldType get ( const int localRow,
                                        const int localCol ) const
      {
        return localMatrix().get( localRow, localCol );
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::scale */
      void scale ( const RangeFieldType& scalar )
      {
        return localMatrix().scale( scalar );
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::clear */
      void clear ()
      {
        return localMatrix().clear();
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::resort */
      void resort ()
      {
        return localMatrix().resort();
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::rows */
      int rows () const
      {
        return localMatrix().rows();
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::columns */
      int columns () const
      {
        return localMatrix().columns();
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::multiplyAdd */
      template <class DomainLocalFunctionImp,
                class RangeLocalFunctionImp>
      void multiplyAdd(const DomainLocalFunctionImp& dLf,
                              RangeLocalFunctionImp& rLf)
      {
        localMatrix().multiplyAdd( dLf, rLf);
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::domainSpace */
      const DomainSpaceType &domainSpace () const
      {
        return localMatrix().domainSpace();
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::rangeSpace */
      const RangeSpaceType &rangeSpace () const
      {
        return localMatrix().rangeSpace();
      }

      const DomainEntityType &domainEntity () const { return localMatrix().domainEntity(); }
      const RangeEntityType &rangeEntity () const { return localMatrix().rangeEntity(); }

      /** \copydoc Dune::Fem::LocalMatrixInterface::domainBasisFunctionSet */
      const DomainBasisFunctionSetType &domainBasisFunctionSet () const
      {
        return localMatrix().domainBasisFunctionSet();
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::rangeBasisFunctionSet */
      const RangeBasisFunctionSetType &rangeBasisFunctionSet () const
      {
        return localMatrix().rangeBasisFunctionSet();
      }

    protected:
      const WrappedLocalMatrixType &localMatrix () const { return localMatrix_; }
      WrappedLocalMatrixType &localMatrix () { return localMatrix_; }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_LOCALMATRIXWRAPPER_HH
