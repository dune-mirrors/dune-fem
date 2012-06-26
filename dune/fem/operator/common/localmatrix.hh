#ifndef DUNE_LOCALMATRIX_HH
#define DUNE_LOCALMATRIX_HH

//- Dune includes 
#include <dune/fem/misc/bartonnackmaninterface.hh>

namespace Dune 
{ 

  namespace Fem
  {

  /** @ingroup Matrix  
      @{ 
  **/

  //- forward declaration of MatrixColumnObject, see below 
  template <class Traits> 
  class MatrixColumnObject ;

  /** \brief Interface for local matrix classes. */
  template< class LocalMatrixTraits >
  class LocalMatrixInterface 
  : public BartonNackmanInterface< LocalMatrixInterface< LocalMatrixTraits >,
                                   typename LocalMatrixTraits::LocalMatrixType >
  {
    typedef LocalMatrixInterface< LocalMatrixTraits > ThisType;
    typedef BartonNackmanInterface< LocalMatrixInterface< LocalMatrixTraits >,
                                    typename LocalMatrixTraits::LocalMatrixType >
      BaseType;

  public:  
    //! type of traits class 
    typedef LocalMatrixTraits Traits; 

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

    /*! type of block (i.e. FieldMatrix for BlockMatrices */
    typedef typename Traits :: LittleBlockType  LittleBlockType;

    typedef MatrixColumnObject< Traits >  MatrixColumnType ;
  protected:  
    using BaseType::asImp;

    //! constructor 
    LocalMatrixInterface ()
    {}
   
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
    inline void add ( const int localRow, 
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

    /*! \brief set row to unit rwo (all entries zero only diagonal entry is one)
        \param[in] localRow local row that is set to unit row 
    */
    /* makes no sens for unsymetrix matrices, use clearRow( localRow ) and 
     * set( localRow, localRow, 1.) instead */
    inline void unitRow ( const int localRow ) DUNE_DEPRECATED
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().unitRow( localRow ));
    }


    /*! \brief set row to zero values
        \param[in] localRow local row that is set to zero 
     */
    inline void clearRow( const int  localRow )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
          asImp().clearRow( localRow ));
    }

    /*! \brief ser column entries to zero
        \param[in] localCol local column that is set to zero
     */

    inline void clearCol( const int localCol )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
          asImp().clearCol( localCol ));
    }   
    

    /*! \brief multiply left hand side with local matrix and add to right hand side 
               rhs += Matrix * lhs  
        \param[in] lhs left hand side 
        \param[out] rhs right hand side 
    */
    template <class DomainLocalFunctionType,
              class RangeLocalFunctionType> 
    inline void multiplyAdd(const DomainLocalFunctionType& lhs,
                            RangeLocalFunctionType& rhs) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().multiplyAdd( lhs, rhs ) );
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

    /*! \brief scale matrix with scalar value 
        \param[in] scalar scalar value that scales the matrix
    */
    inline void scale ( const RangeFieldType& scalar )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().scale( scalar ) );
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
    inline const RangeBaseFunctionSetType &rangeBaseFunctionSet () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().rangeBaseFunctionSet() );
      return asImp().rangeBaseFunctionSet();
    }

    /** \brief return column object for local matrix which contains axpy methods 
               for convenience 
        \param col  local column number 

        \return object of type MatrixColumnObject 
    */
    MatrixColumnType column( const unsigned int col )  
    {
      return MatrixColumnType( asImp(), col );
    }
  };


  /** \brief Default implementation for local matrix classes. */
  template< class LocalMatrixTraits >
  class LocalMatrixDefault
  : public LocalMatrixInterface< LocalMatrixTraits >
  {
    typedef LocalMatrixDefault< LocalMatrixTraits > ThisType;
    typedef LocalMatrixInterface< LocalMatrixTraits > BaseType;

  public:
    typedef LocalMatrixTraits Traits;

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
    LocalMatrixDefault ( const DomainSpaceType &domainSpace,
                         const RangeSpaceType &rangeSpace )
    : domainSpace_( domainSpace ),
      rangeSpace_( rangeSpace ),
      domainBaseSet_(),
      rangeBaseSet_()
    {}

    template< class DomainEntityType, class RangeEntityType >
    LocalMatrixDefault ( const DomainSpaceType &domainSpace,
                         const RangeSpaceType &rangeSpace,
                         const DomainEntityType &domainEntity,
                         const RangeEntityType &rangeEntity )
    : domainSpace_( domainSpace ),
      rangeSpace_( rangeSpace ),
      domainBaseSet_( domainSpace.baseFunctionSet( domainEntity ) ),
      rangeBaseSet_( rangeSpace.baseFunctionSet( rangeEntity ) )
    {}

    LocalMatrixDefault ( const LocalMatrixDefault& org )
    : domainSpace_( org.domainSpace_ ),
      rangeSpace_( org.rangeSpace_ ),
      domainBaseSet_( org.domainBaseSet_ ),
      rangeBaseSet_( org.rangeBaseSet_ )
    {}

  public:
    /** \copydoc Dune::LocalMatrixInterface::init */
    template< class DomainEntityType, class RangeEntityType >
    inline void init ( const DomainEntityType &domainEntity,
                       const RangeEntityType &rangeEntity )
    {
      domainBaseSet_ = domainSpace_.baseFunctionSet( domainEntity );
      rangeBaseSet_ = rangeSpace_.baseFunctionSet( rangeEntity );
    }

    /** \copydoc Dune::LocalMatrixInterface::resort */
    void resort () {}

    /** \copydoc Dune::LocalMatrixInterface::rows */
    int rows () const { return rangeBaseSet_.size(); }
    
    /** \copydoc Dune::LocalMatrixInterface::columns */
    int columns () const { return domainBaseSet_.size(); }

    /** \copydoc Dune::LocalMatrixInterface::domainSpace */
    const DomainSpaceType &domainSpace () const { return domainSpace_; }
    
    /** \copydoc Dune::LocalMatrixInterface::rangeSpace */
    const RangeSpaceType &rangeSpace () const { return rangeSpace_; }

    /** \copydoc Dune::LocalMatrixInterface::domainBaseFunctionSet */
    const DomainBaseFunctionSetType &domainBaseFunctionSet () const
    {
      return domainBaseSet_;
    }
    
    /** \copydoc Dune::LocalMatrixInterface::rangeBaseFunctionSet */
    const RangeBaseFunctionSetType &rangeBaseFunctionSet () const
    {
      return rangeBaseSet_;
    }

    /** \copydoc Dune::LocalMatrixInterface::multiplyAdd */
    template <class DomainLocalFunctionType,
              class RangeLocalFunctionType> 
    inline void multiplyAdd(const DomainLocalFunctionType& lhs,
                            RangeLocalFunctionType& rhs) const
    {
      const int row = this->rows();
      const int col = this->columns();
      for(int i=0; i<row; ++i)
      {
        for(int j=0; j<col; ++j)
        {
          rhs[i] += this->get(i,j) * lhs[j];
        }
      }
    }
  };

  template <class Traits> 
  class MatrixColumnObject 
  {
  public:
    //! type of local matrix implementation 
    typedef typename Traits :: LocalMatrixType LocalMatrixType;

    //! type of domain discrete function space
    typedef typename Traits :: DomainSpaceType  DomainSpaceType;

    //! type of range 
    typedef typename DomainSpaceType :: RangeType          RangeType ;
    //! type of jacobian range 
    typedef typename DomainSpaceType :: JacobianRangeType  JacobianRangeType ;
    //! type of range field 
    typedef typename DomainSpaceType :: RangeFieldType     RangeFieldType ;

  protected:
    // reference to local matrix 
    LocalMatrixType& localMatrix_;
    // local column number 
    const unsigned int column_;

    //! constructor taking local matrix and column number 
    MatrixColumnObject( LocalMatrixType& localMatrix, const unsigned int col ) 
      : localMatrix_( localMatrix ),
        column_( col ) 
    {
    }

    // at the moment only allow LocalMatrixInterface to construct this object 
    friend class LocalMatrixInterface< Traits >;

  public:

    /** \brief axpy operation for local matrices 
     *
     *  Denoting an entry of the local matrix by \f$a_{i,j}\f$ and the base
     *  functions by \f$\varphi_i\f$, this function performs the following
     *  operation:
     *  \f[
     *  a_{i,j} = a_{i,j} + weight * (factor \cdot \varphi_i( x ))
     *  \f]
     *  \param[in]  phi     evaluations of all base functions \f$\varphi_i( x )\f$
     *  \param[in]  factor  axpy factor
     *  \param[in]  weight  integration weight for quadrature point (default = 1)
     */
    template <class RangeVectorType>
    void axpy( const RangeVectorType& phi, 
               const RangeType& factor,
               const RangeFieldType& weight = RangeFieldType(1) ) 
    {
      const unsigned int numBaseFunctions = localMatrix_.rows();
      assert( phi.size() >= numBaseFunctions );
      for( unsigned int row = 0; row < numBaseFunctions; ++ row )
      {
        RangeFieldType value = factor * phi[ row ];
        localMatrix_.add( row, column_, weight * value );
      }
    }

    /** \brief axpy operation for local matrices 
     *
     *  Denoting an entry of the local matrix by \f$a_{i,j}\f$ and 
     *  the gradients of the  base functions by \f$\nabla \varphi_i\f$, 
     *  this function performs the following operation:
     *  \f[
     *  a_{i,j} = a_{i,j} + weight * (jacobianFactor \cdot \nabla\varphi_i( x ))
     *  \f]
     *  \param[in]  dphi           evaluations of the jacobian of all base functions \f$\varphi_i( x )\f$
     *  \param[in]  jacobianFactor axpy factor
     *  \param[in]  weight         integration weight for quadrature point (default = 1)
     */
    template <class JacobianVectorType>
    void axpy( const JacobianVectorType& dphi,
               const JacobianRangeType& jacobianFactor,
               const RangeFieldType& weight = RangeFieldType(1) ) 
    {
      const unsigned int numBaseFunctions = localMatrix_.rows();
      assert( dphi.size() >= numBaseFunctions );
      for( unsigned int row = 0; row < numBaseFunctions; ++ row )
      {
        RangeFieldType value = 0; 
        for( int k = 0; k < jacobianFactor.rows; ++k )
          value += jacobianFactor[ k ] * dphi[ row ][ k ];

        localMatrix_.add( row, column_, weight * value );
      }
    }

    /** \brief axpy operation for local matrices 
     *
     *  Denoting an entry of the local matrix by \f$a_{i,j}\f$ and 
     *  the base functions by \f$\nabla \varphi_i\f$, 
     *  this function performs the following operation:
     *  \f[
     *  a_{i,j} = a_{i,j} + weight (factor \cdot \varphi_i( x ) + jacobianFactor \cdot \nabla\varphi_i( x ))
     *  \f]
     *  \param[in]  phi            evaluations of all base functions \f$\varphi_i( x )\f$
     *  \param[in]  dphi           evaluations of the jacobian of all base functions \f$\varphi_i( x )\f$
     *  \param[in]  factor         axpy factor for phi 
     *  \param[in]  jacobianFactor axpy factor for dphi 
     *  \param[in]  weight         integration weight for quadrature point (default = 1)
     */
    template <class RangeVectorType, class JacobianVectorType>
    void axpy( const RangeVectorType& phi, 
               const JacobianVectorType& dphi,
               const RangeType& factor,
               const JacobianRangeType& jacobianFactor,
               const RangeFieldType& weight = RangeFieldType(1) ) 
    {
      const unsigned int numBaseFunctions = localMatrix_.rows();
      assert( phi.size() >= numBaseFunctions );
      assert( dphi.size() >= numBaseFunctions );
      for( unsigned int row = 0; row < numBaseFunctions; ++ row )
      {
        RangeFieldType value = factor * phi[ row ];
        for( int k = 0; k < jacobianFactor.rows; ++k )
          value += jacobianFactor[ k ] * dphi[ row ][ k ];

        localMatrix_.add( row, column_, weight * value );
      }
    }
  };

///@}

  } // end namespace Fem 

  // #if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: LocalMatrixDefault ;
  using Fem :: MatrixColumnObject ;

  // #endif // DUNE_FEM_COMPATIBILITY

} // end namespace Dune 

#endif // #ifndef DUNE_LOCALMATRIX_HH
