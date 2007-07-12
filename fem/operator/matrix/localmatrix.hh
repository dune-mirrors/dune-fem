#ifndef DUNE_FEM_LOCALMATRIX_HH
#define DUNE_FEM_LOCALMATRIX_HH

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

  template< class TraitsImp, class LocalMatrixImp >
  class LocalMatrixInterface
  {
  public:
    typedef TraitsImp TraitsType;
    typedef typename TraitsType :: RowType RowType;
    typedef typename TraitsType :: DomainFunctionSpaceType DomainFunctionSpaceType;
    typedef typename TraitsType :: RangeFunctionSpaceType RangeFunctionSpaceType;

    typedef typename RowType :: FieldType FieldType;
    typedef typename DomainFunctionSpaceType :: BaseFunctionSetType
      DomainBaseFunctionSetType;
    typedef typename RangeFunctionSpaceType :: BaseFunctionSetType
      RangeBaseFunctionSetType;

  private:
    typedef LocalMatrixInterface< TraitsType, LocalMatrixImp > ThisType;

  public:
    typedef ThisType LocalMatrixInterfaceType;

  public:
    inline LocalMatrixImp &operator= ( const FieldType &s )
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().operator=( s ) );
      return asImp();
    }

    inline const RowType operator[] ( unsigned int row ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp()[ row ] );
      return asImp()[ row ];
    }

    inline RowType operator[] ( unsigned int row )
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp()[ row ] );
      return asImp()[ row ];
    }

    inline unsigned int columns () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().columns() );
      return asImp().columns();
    }

    inline const DomainBaseFunctionSetType &domainBaseFunctionSet () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().domainBaseFunctionSet() );
      return asImp().domainBaseFunctionSet();
    }

    inline const RangeBaseFunctionSetType &rangeBaseFunctionSet () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().rangeBaseFunctionSet() );
      return asImp().rangeBaseFunctionSet();
    }

    inline unsigned int rows () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().rows() );
      return asImp().rows();
    }

  protected:
    inline const LocalMatrixImp &asImp () const
    {
      return static_cast< const LocalMatrixImp& >( *this );
    }

    inline LocalMatrixImp &asImp ()
    {
      return static_cast< LocalMatrixImp& >( *this );
    }
  };
  
 
  
  template< class TraitsImp, class LocalMatrixImp >
  class LocalMatrixDefault
  : public LocalMatrixInterface< TraitsImp, LocalMatrixImp >
  {
  public:
    typedef TraitsImp TraitsType;
    typedef typename TraitsType :: RowType RowType;
    typedef typename TraitsType :: DomainFunctionSpaceType DomainFunctionSpaceType;
    typedef typename TraitsType :: RangeFunctionSpaceType RangeFunctionSpaceType;

    typedef typename RowType :: FieldType FieldType;
    typedef typename DomainFunctionSpaceType :: BaseFunctionSetType
      DomainBaseFunctionSetType;
    typedef typename RangeFunctionSpaceType :: BaseFunctionSetType
      RangeBaseFunctionSetType;

  private:
    typedef LocalMatrixDefault< TraitsType, LocalMatrixImp > ThisType;
    typedef LocalMatrixInterface< TraitsType, LocalMatrixImp > BaseType;

    using BaseType :: asImp;

  public:
    typedef ThisType LocalMatrixInterfaceType;

  protected:
    const DomainFunctionSpaceType &domainFunctionSpace_;
    const RangeFunctionSpaceType &rangeFunctionSpace_;

    DomainBaseFunctionSetType domainBaseFunctionSet_;
    RangeBaseFunctionSetType rangeBaseFunctionSet_;
    
  public:
    inline LocalMatrixDefault ( const DomainFunctionSpaceType &domainFunctionSpace,
                                const RangeFunctionSpaceType &rangeFunctionSpace )
    : domainFunctionSpace_( domainFunctionSpace ),
      rangeFunctionSpace_( rangeFunctionSpace ),
      domainBaseFunctionSet_(),
      rangeBaseFunctionSet_()
    {
    }

    template< class EntityType >
    inline LocalMatrixDefault ( const DomainFunctionSpaceType &domainFunctionSpace,
                                const RangeFunctionSpaceType &rangeFunctionSpace,
                                const EntityType &entity )
    : domainFunctionSpace_( domainFunctionSpace ),
      rangeFunctionSpace_( rangeFunctionSpace ),
      domainBaseFunctionSet_( domainFunctionSpace_.baseFunctionSet( entity ) ),
      rangeBaseFunctionSet_( rangeFunctionSpace_.baseFunctionSet( entity ) )
    {
    }

    inline LocalMatrixImp &operator= ( const FieldType &s )
    {
      const unsigned int rows = asImp().rows();
      for( unsigned int i = 0; i < rows; ++i )
        (*this)[ i ] = s;
      return asImp();
    }

    inline unsigned int columns () const
    {
      return domainBaseFunctionSet_.numBaseFunctions();
    }

    inline const DomainBaseFunctionSetType &domainBaseFunctionSet () const
    {
      return domainBaseFunctionSet_;
    }

    template< class EntityType >
    inline void init ( const EntityType &entity )
    {
      domainBaseFunctionSet_ = domainFunctionSpace_.baseFunctionSet( entity );
      rangeBaseFunctionSet_ = rangeFunctionSpace_.baseFunctionSet( entity );
    }

    inline const RangeBaseFunctionSetType &rangeBaseFunctionSet () const
    {
      return rangeBaseFunctionSet_;
    }

    inline unsigned int rows () const
    {
      return rangeBaseFunctionSet_.numBaseFunctions();
    }
  };



  template< class DomainFunctionSpaceImp, class RangeFunctionSpaceImp >
  class TemporaryLocalMatrixRow;



  template< class DomainFunctionSpaceImp, class RangeFunctionSpaceImp >
  class TemporaryLocalMatrixTraits
  {
  public:
    typedef DomainFunctionSpaceImp DomainFunctionSpaceType;
    typedef RangeFunctionSpaceImp RangeFunctionSpaceType;

  private:
    typedef TemporaryLocalMatrixTraits< DomainFunctionSpaceType, RangeFunctionSpaceType >
      ThisType;

  public:
    typedef TemporaryLocalMatrixRow< DomainFunctionSpaceType, RangeFunctionSpaceType >
      RowType;
  };

  

  template< class DomainFunctionSpaceImp, class RangeFunctionSpaceImp >
  class TemporaryLocalMatrix
  : public LocalMatrixDefault< TemporaryLocalMatrixTraits< DomainFunctionSpaceImp,
                                                           RangeFunctionSpaceImp >,
                               TemporaryLocalMatrix< DomainFunctionSpaceImp,
                                                     RangeFunctionSpaceImp > >
  {
  public:
    typedef DomainFunctionSpaceImp DomainFunctionSpaceType;
    typedef RangeFunctionSpaceImp RangeFunctionSpaceType;

    typedef TemporaryLocalMatrixTraits< DomainFunctionSpaceType, RangeFunctionSpaceType >
      TraitsType;

    typedef typename TraitsType :: RowType RowType;
    typedef typename RowType :: FieldType FieldType;

  private:
    typedef TemporaryLocalMatrix< DomainFunctionSpaceType, RangeFunctionSpaceType >
      ThisType;
    typedef LocalMatrixDefault< TraitsType, ThisType > BaseType;

    friend class TemporaryLocalMatrixRow< DomainFunctionSpaceType, RangeFunctionSpaceType >;

    using BaseType :: domainBaseFunctionSet_;
    using BaseType :: rangeBaseFunctionSet_;

  protected:
    unsigned int columns_;
    unsigned int rows_;

    unsigned int size_;
    FieldType *fields_;

  public:
    inline TemporaryLocalMatrix ( const DomainFunctionSpaceType &domainFunctionSpace,
                                  const RangeFunctionSpaceType &rangeFunctionSpace )
    : BaseType( domainFunctionSpace, rangeFunctionSpace ),
      columns_( 0 ),
      rows_( 0 ),
      size_( 0 ),
      fields_( NULL )
    {
    }
    
    template< class EntityType >
    inline TemporaryLocalMatrix ( const DomainFunctionSpaceType &domainFunctionSpace,
                                  const RangeFunctionSpaceType &rangeFunctionSpace,
                                  const EntityType &entity )
    : BaseType( domainFunctionSpace, rangeFunctionSpace, entity ),
      columns_( domainBaseFunctionSet_.numBaseFunctions() ),
      rows_( rangeBaseFunctionSet_.numBaseFunctions() ),
      size_( rows_ * columns_ ),
      fields_( new FieldType[ size_ ] )
    {
      assert( fields_ != NULL );
    }

    inline ~TemporaryLocalMatrix ()
    {
      if( fields_ != NULL )
        delete[]( fields_ );
    }

    inline ThisType &operator= ( const FieldType &s )
    {
      return BaseType :: operator=( s );
    }

    inline const RowType operator[] ( unsigned int row ) const
    {
      return RowType( *this, row );
    }
    
    inline RowType operator[] ( unsigned int row )
    {
      return RowType( *this, row );
    }

    inline unsigned int columns () const
    {
      return columns_;
    }

    template< class EntityType >
    inline void init ( const EntityType &entity )
    {
      BaseType :: init( entity );

      columns_ = domainBaseFunctionSet_.numBaseFunctions();
      rows_ = rangeBaseFunctionSet_.numBaseFunctions();
      const unsigned int size = rows_ * columns_;

      if( size > size_ )
      {
        size_ = size;
        if( fields_ != NULL )
          delete[] fields_;
        fields_ = new FieldType[ size_ ];
        
        assert( fields_ != NULL );
      }
    }
    
    inline unsigned int rows () const
    {
      return rows_;
    }
  };


  
  template< class DomainFunctionSpaceImp, class RangeFunctionSpaceImp >
  class TemporaryLocalMatrixRow
  {
  public:
    typedef DomainFunctionSpaceImp DomainFunctionSpaceType;
    typedef RangeFunctionSpaceImp RangeFunctionSpaceType;

  private:
    typedef TemporaryLocalMatrixRow< DomainFunctionSpaceType, RangeFunctionSpaceType >
      ThisType;

    typedef TemporaryLocalMatrix< DomainFunctionSpaceType, RangeFunctionSpaceType >
      TemporaryLocalMatrixType;

  public:
    typedef typename RangeFunctionSpaceType :: RangeFieldType FieldType;

  protected:
    TemporaryLocalMatrixType &matrix_;
    unsigned int index_;
    
  public:
    inline TemporaryLocalMatrixRow ( TemporaryLocalMatrixType &matrix,
                                     unsigned int row )
    : matrix_( matrix ),
      index_( row * matrix_.columns_ )
    {
      assert( row < matrix_.rows_ );
    }

    inline ThisType &operator= ( const FieldType &s )
    {
      for( unsigned int i = 0; i < matrix_.columns_; ++i )
        matrix_.fields_[ index_ + i ] = s;
      return *this;
    }

    inline const FieldType &operator[] ( unsigned int column ) const
    {
      assert( column < matrix_.columns_ );
      return matrix_.fields_[ index_ + column ];
    }

    inline FieldType &operator[] ( unsigned int column )
    {
      assert( column < matrix_.columns_ );
      return matrix_.fields_[ index_ + column ];
    }

    inline unsigned int size () const
    {
      return matrix_.columns_;
    }
  };

}

#endif
