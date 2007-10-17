#ifndef DUNE_FEM_LOCALMATRIXWRAPPER_HH
#define DUNE_FEM_LOCALMATRIXWRAPPER_HH

#include <dune/fem/operator/common/localmatrix.hh>

namespace Dune
{

  template< class LocalMatrixStackImp >
  class LocalMatrixWrapper;



  template< class LocalMatrixStackImp >
  struct LocalMatrixWrapperTraits
  {
    typedef LocalMatrixStackImp LocalMatrixStackType;

    typedef typename LocalMatrixStackImp :: ObjectType WrappedLocalMatrixType;

    typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;

    typedef typename WrappedLocalMatrixType :: RangeFieldType RangeFieldType;
    
    typedef typename WrappedLocalMatrixType :: DomainSpaceType DomainSpaceType;
    typedef typename WrappedLocalMatrixType :: RangeSpaceType RangeSpaceType;

    typedef typename WrappedLocalMatrixType :: LittleBlockType LittleBlockType;
  };



  template< class LocalMatrixStackImp >
  class LocalMatrixWrapper
  : LocalMatrixInterface< LocalMatrixWrapperTraits< LocalMatrixStackImp > >
  {
  public:
    //! type of the local matrix stack
    typedef LocalMatrixStackImp LocalMatrixStackType;

    //! type of the traits
    typedef LocalMatrixWrapperTraits< LocalMatrixStackType > Traits;

  private:
    typedef LocalMatrixWrapper< LocalMatrixStackType > ThisType;
    typedef LocalMatrixInterface< Traits > BaseType;

  public:
    //! type of the wrapped local matrix
    typedef typename Traits :: WrappedLocalMatrixType WrappedLocalMatrixType;

    typedef typename Traits :: RangeFieldType RangeFieldType;

    typedef typename BaseType :: DomainSpaceType DomainSpaceType;
    typedef typename BaseType :: RangeSpaceType RangeSpaceType;

    typedef typename BaseType :: DomainBaseFunctionSetType
      DomainBaseFunctionSetType;
    typedef typename BaseType :: RangeBaseFunctionSetType
      RangeBaseFunctionSetType;

  private:
    typedef typename LocalMatrixStackType :: PointerType
      WrappedLocalMatrixPtrType;

  private:
    // ObjectPointer to the actual local matrix
    // (the pointer is required to keep the reference alive)
    WrappedLocalMatrixPtrType localMatrixPtr_;

    // reference to the actual local matrix
    WrappedLocalMatrixType &localMatrix_;

  public:
    //! constructor creating an uninitialized local matrix
    inline explicit LocalMatrixWrapper ( LocalMatrixStackType &stack )
    : localMatrixPtr_( stack.getObject() ),
      localMatrix_( *localMatrixPtr_ )
    {
    }
    
    //! constructor initializing the wrapped local matrix
    template< class DomainEntityType, class RangeEntityType >
    inline LocalMatrixWrapper( LocalMatrixStackType &stack,
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
    inline LocalMatrixWrapper ( const ThisType &other )
    : localMatrixPtr_( other.localMatrixPtr_ ),
      localMatrix_( *localMatrixPtr_ )
    {
    }

  private:
    // prohibit assignment
    inline ThisType &operator= ( const ThisType & );

  public:
    /** \copydoc Dune::LocalMatrixInterface::init */
    template< class DomainEntityType, class RangeEntityType >
    inline void init ( const DomainEntityType &domainEntity,
                       const RangeEntityType &rangeEntity )
    {
      localMatrix().init( domainEntity, rangeEntity );
    }

    /** \copydoc Dune::LocalMatrixInterface::add */
    inline void add ( const int localRow,
                      const int localCol,
                      const RangeFieldType &value )
    {
      localMatrix().add( localRow, localCol, value );
    }
    
    /** \copydoc Dune::LocalMatrixInterface::set */
    inline void set ( const int localRow,
                      const int localCol,
                      const RangeFieldType &value )
    {
      localMatrix().set( localRow, localCol, value );
    }
    
    /** \copydoc Dune::LocalMatrixInterface::unitRow */
    inline void unitRow ( const int localRow )
    {
      localMatrix().unitRow( localRow ); 
    }
    
    /** \copydoc Dune::LocalMatrixInterface::get */
    inline const RangeFieldType get ( const int localRow,
                                      const int localCol ) const
    {
      return localMatrix().get( localRow, localCol );
    }

    /** \copydoc Dune::LocalMatrixInterface::clear */
    inline void clear ()
    {
      return localMatrix().clear();
    }

    /** \copydoc Dune::LocalMatrixInterface::resort */
    inline void resort ()
    {
      return localMatrix().resort();
    }

    /** \copydoc Dune::LocalMatrixInterface::rows */
    inline int rows () const
    {
      return localMatrix().rows();
    }

    /** \copydoc Dune::LocalMatrixInterface::columns */
    inline int columns () const
    {
      return localMatrix().columns();
    }

    /** \copydoc Dune::LocalMatrixInterface::multiplyAdd */
    template <class DomainLocalFunctionImp, 
              class RangeLocalFunctionImp>
    inline void multiplyAdd(const DomainLocalFunctionImp& dLf,
                            RangeLocalFunctionImp& rLf)
    {
      localMatrix().multiplyAdd( dLf, rLf); 
    }
    
    /** \copydoc Dune::LocalMatrixInterface::domainSpace */
    inline const DomainSpaceType &domainSpace () const
    {
      return localMatrix().domainSpace();
    }

    /** \copydoc Dune::LocalMatrixInterface::rangeSpace */
    inline const RangeSpaceType &rangeSpace () const
    {
      return localMatrix().rangeSpace();
    }
    
    /** \copydoc Dune::LocalMatrixInterface::domainBaseFunctionSet */
    inline const DomainBaseFunctionSetType &domainBaseFunctionSet () const
    {
      return localMatrix().domainBaseFunctionSet();
    }

    /** \copydoc Dune::LocalMatrixInterface::rangeBaseFunctionSet */
    inline const RangeBaseFunctionSetType &rangeBaseFunctionSet () const
    {
      return localMatrix().rangeBaseFunctionSet();
    }

  protected:
    inline const WrappedLocalMatrixType &localMatrix () const
    {
      return localMatrix_;
    }

    inline WrappedLocalMatrixType &localMatrix ()
    {
      return localMatrix_;
    }
  };
  
}

#endif
