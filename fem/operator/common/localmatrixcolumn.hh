#ifndef DUNE_FEM_OPERATOR_COMMON_LOCALMATRIXCOLUMN_HH
#define DUNE_FEM_OPERATOR_COMMON_LOCALMATRIXCOLUMN_HH

#include <utility>

namespace Dune
{

  namespace Fem
  {

    // LocalMatrixEntry
    // ----------------

    template< class LocalMatrix >
    class LocalMatrixEntry
    {
      typedef LocalMatrixEntry< LocalMatrix > ThisType;

    public:
      typedef LocalMatrix LocalMatrixType;

      typedef typename LocalMatrixType::RangeFieldType RangeFieldType;

      LocalMatrixEntry ( const LocalMatrixType &localMatrix, unsigned int row, unsigned int column )
        : localMatrix_( localMatrix ), row_( row ), column_( column )
      {}

      operator RangeFieldType () const { return localMatrix_.get( row_, column_ ); }

      ThisType &operator= ( const RangeFieldType &value ) { localMatrix_.set( row_, col_, value ); return *this; }

      ThisType &operator+= ( const RangeFieldType &value ) { localMatrix_.add( row_, col_, value ); return *this; }
      ThisType &operator-= ( const RangeFieldType &value ) { localMatrix_.add( row_, col_, -value ); return *this; }

      ThisType &operator*= ( const RangeField &value ) { localMatrix_.set( row_, col_, localMatrix_.get( row_, col_ ) * value ); return *this; }
      ThisType &operator/= ( const RangeField &value ) { localMatrix_.set( row_, col_, localMatrix_.get( row_, col_ ) / value ); return *this; }

    private:
      LocalMatrixType &localMatrix_;
      unsigned int row_, column_;
    };



    // LocalMatrixColumn
    // -----------------

    template< class LocalMatrix >
    class LocalMatrixColumn
    {
      typedef LocalMatrixColumn< LocalMatrix > ThisType;

    public:
      typedef LocalMatrix LocalMatrixType;

      typedef typename LocalMatrixType::RangeFieldType RangeFieldType;

      typedef typename LocalMatrixType::RangeType RangeType;
      typedef typename LocalMatrixType::JacobianRangeType JacobianRangeType;

      typedef typename LocalMatrixType::RangeBasisFunctionSetType BasisFunctionSetType;

      LocalMatrixColumn ( LocalMatrixType &localMatrix, unsigned int column )
        : localMatrix_( localMatrix ), column_( column )
      {}

      RangeFieldType operator[] ( unsigned int row ) const { return localMatrix_.get( row, column_ ); }

      LocalMatrixEntry< LocalMatrix > operator[] ( unsigned int row ) { return LocalMatrixEntry( localMatrix_, row, col_ ); }

      template< class Point, class... Factor >
      void axpy ( const Point &x, Factor &&... factor )
      {
        basisFunctionSet().axpy( x, std::forward< Factor >( factor )..., *this );
      }

      BasisFunctionSetType &basisFunctionSet () const { return localMatrix_.rangeBasisFunctionSet(); }

    private:
      LocalMatrixType &localMatrix_;
      unsigned int column_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_COMMON_LOCALMATRIXCOLUMN_HH
