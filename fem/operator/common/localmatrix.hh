#ifndef DUNE_LOCALMATRIX_HH
#define DUNE_LOCALMATRIX_HH

//- system includes 
//#include <vector> -- whatfor? 

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune 
{ 

  /** @ingroup Matrix  
      @{ 
  **/

  /** \brief Interface for local matrix classes. */
  template< class LocalMatrixTraits >
  class LocalMatrixInterface 
  {
  public:  
    //! type of traits class 
    typedef LocalMatrixTraits Traits; 
   
  private:
    typedef LocalMatrixInterface< Traits > ThisType;

  public:
    //! type of this interface
    typedef ThisType LocalMatrixInterfaceType;

    //! type of local matrix implementation 
    typedef typename Traits :: LocalMatrixType LocalMatrixType; 

    //! type of range field 
    typedef typename Traits :: RangeFieldType RangeFieldType;

    //! type of domain discrete function space
    typedef typename Traits :: DomainSpaceType DomainSpaceType;
    
    //! type of range discrete function space
    typedef typename Traits :: RangeSpaceType RangeSpaceType;

    //! type of base function sets within domain function space
    typedef typename DomainSpaceType :: BaseFunctionSetType
      DomainBaseFunctionSetType;

    //! type of base function sets within range function space
    typedef typename RangeSpaceType :: BaseFunctionSetType
      RangeBaseFunctionSetType;

    /* Is such a declaration valid?
    //! 
    typedef typename LocalMatrixType:: block_type LittleBlockType;
    */

  protected:  
    //! constructor 
    inline LocalMatrixInterface ()
    {
    }
   
  public:
    /** \brief initialize the local matrix to entities
     *  \param[in]  domainEntity  entity within grid of domain space,
     *  \param[in]  rangeEntity   entity within grid of range space
     */
    template< class DomainEntityType, class RangeEntityType >
    inline void init ( const DomainEntityType &domainEntity,
                       const RangeEntityType &rangeEntity )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().init( domainEntity, rangeEntity ) );
    }

    /*! \brief add value to matrix entry (row,col) where row and col are
        local row and local column 
        \param[in] localRow local row 
        \param[in] localCol local column 
        \param[in] value value to add 
    */
    inline void add( const int localRow, 
                     const int localCol, 
                     const RangeFieldType &value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().add(localRow,localCol,value));
    }

    /*! \brief set value of matrix entry (row,col) where row and col are
        local row and local column 
        \param[in] localRow local row 
        \param[in] localCol local column 
        \param[in] value value to set  
    */
    inline void set ( const int localRow, 
                      const int localCol, 
                      const RangeFieldType &value )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().set(localRow,localCol,value));
    }

    /*! \brief get value of matrix entry (row,col) where row and col are
        local row and local column 
        \param[in] localRow local row 
        \param[in] localCol local column 
        \return value of matrix entry   
    */
    inline const RangeFieldType get ( const int localRow, 
                                      const int localCol ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(
        asImp().get(localRow,localCol));
      return asImp().get(localRow,localCol);
    }

    /*! \brief set all entries of local matrix to zero */
    inline void clear ()
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().clear());
    }

    /*! \brief resort ordering in global matrix (if possible) */
    inline void resort ()
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().resort());
    }

    /** \brief get number of rows within the matrix */
    inline int rows () const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().rows() );
      return asImp().rows();
    }
    
    /** \brief get number of columns within the matrix */
    inline int columns () const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().columns() );
      return asImp().columns();
    }

    /** \brief access to the domain space */
    inline const DomainSpaceType &domainSpace () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().domainSpace() );
      return asImp().domainSpace();
    }
    
    /** \brief access to the range space */
    inline const RangeSpaceType &rangeSpace () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().rangeSpace() );
      return asImp().rangeSpace();
    }

    /** \brief access to the base function set within the domain space */
    inline const DomainBaseFunctionSetType &domainBaseFunctionSet () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().domainBaseFunctionSet() );
      return asImp().domainBaseFunctionSet();
    }
    
    /** \brief access to the base function set within the range space */
    inline const DomainBaseFunctionSetType &rangeBaseFunctionSet () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().rangeBaseFunctionSet() );
      return asImp().rangeBaseFunctionSet();
    }
    
  private:
    //! Barton-Nackman Trick
    inline LocalMatrixType &asImp ()
    {
      return static_cast<LocalMatrixType&> (*this);
    }
    
    //! Barton-Nackman Trick
    inline const LocalMatrixType &asImp () const
    {
      return static_cast<const LocalMatrixType&> (*this);
    }
  };
 


  /** \brief Default implementation for local matrix classes. */
  template< class LocalMatrixTraits >
  class LocalMatrixDefault
  : public LocalMatrixInterface< LocalMatrixTraits >
  {
  public:
    typedef LocalMatrixTraits Traits;

  private:
    typedef LocalMatrixDefault< Traits > ThisType;
    typedef LocalMatrixInterface< Traits > BaseType;

  public:
    typedef typename BaseType :: DomainSpaceType DomainSpaceType;
    typedef typename BaseType :: RangeSpaceType RangeSpaceType;

    typedef typename BaseType :: DomainBaseFunctionSetType
      DomainBaseFunctionSetType;
    typedef typename BaseType :: RangeBaseFunctionSetType
      RangeBaseFunctionSetType;

  protected:
    const DomainSpaceType &domainSpace_;
    const RangeSpaceType &rangeSpace_;

    DomainBaseFunctionSetType domainBaseSet_;
    RangeBaseFunctionSetType rangeBaseSet_;

  protected:
    inline LocalMatrixDefault ( const DomainSpaceType &domainSpace,
                                const RangeSpaceType &rangeSpace )
    : domainSpace_( domainSpace ),
      rangeSpace_( rangeSpace ),
      domainBaseSet_(),
      rangeBaseSet_()
    {
    }

    template< class DomainEntityType, class RangeEntityType >
    inline LocalMatrixDefault ( const DomainSpaceType &domainSpace,
                                const RangeSpaceType &rangeSpace,
                                const DomainEntityType &domainEntity,
                                const RangeEntityType &rangeEntity )
    : domainSpace_( domainSpace ),
      rangeSpace_( rangeSpace ),
      domainBaseSet_( domainSpace.baseFunctionSet( domainEntity ) ),
      rangeBaseSet_( rangeSpace.baseFunctionSet( rangeEntity ) )
    {
    }

  public:
    template< class DomainEntityType, class RangeEntityType >
    inline void init ( const DomainEntityType &domainEntity,
                       const RangeEntityType &rangeEntity )
    {
      domainBaseSet_ = domainSpace_.baseFunctionSet( domainEntity );
      rangeBaseSet_ = rangeSpace_.baseFunctionSet( rangeEntity );
    }

    inline void resort ()
    {
    }
    
    inline int rows () const
    {
      return rangeBaseSet_.numBaseFunctions();
    }
    
    inline int columns () const
    {
      return domainBaseSet_.numBaseFunctions();
    }

    inline const DomainSpaceType &domainSpace () const
    {
      return domainSpace_;
    }
    
    inline const RangeSpaceType &rangeSpace () const
    {
      return rangeSpace_;
    }

    inline const DomainBaseFunctionSetType &domainBaseFunctionSet () const
    {
      return domainBaseSet_;
    }
    
    inline const DomainBaseFunctionSetType &rangeBaseFunctionSet () const
    {
      return rangeBaseSet_;
    }
  };

///@} 
} // end namespace Dune 
#endif
