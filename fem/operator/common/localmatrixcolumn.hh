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

      LocalMatrixEntry ( LocalMatrixType &localMatrix, unsigned int row, unsigned int col )
        : localMatrix_( localMatrix ), row_( row ), col_( col )
      {}

      operator RangeFieldType () const { return localMatrix_.get( row_, col_ ); }

      ThisType &operator= ( const RangeFieldType &value ) { localMatrix_.set( row_, col_, value ); return *this; }

      ThisType &operator+= ( const RangeFieldType &value ) { localMatrix_.add( row_, col_, value ); return *this; }
      ThisType &operator-= ( const RangeFieldType &value ) { localMatrix_.add( row_, col_, -value ); return *this; }

      ThisType &operator*= ( const RangeFieldType &value ) { localMatrix_.set( row_, col_, localMatrix_.get( row_, col_ ) * value ); return *this; }
      ThisType &operator/= ( const RangeFieldType &value ) { localMatrix_.set( row_, col_, localMatrix_.get( row_, col_ ) / value ); return *this; }

    private:
      LocalMatrixType &localMatrix_;
      unsigned int row_, col_;
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
      typedef typename LocalMatrixType::RangeBasisFunctionSetType BasisFunctionSetType;

      typedef typename BasisFunctionSetType::RangeType RangeType;
      typedef typename BasisFunctionSetType::JacobianRangeType JacobianRangeType;
      typedef typename BasisFunctionSetType::HessianRangeType HessianRangeType;

      LocalMatrixColumn ( LocalMatrixType &localMatrix, unsigned int col )
        : localMatrix_( localMatrix ), col_( col )
      {}

      RangeFieldType operator[] ( unsigned int row ) const { return localMatrix_.get( row, col_ ); }

      LocalMatrixEntry< LocalMatrixType > operator[] ( unsigned int row ) { return LocalMatrixEntry< LocalMatrixType >( localMatrix_, row, col_ ); }

      template< class Point, class... Factor >
      void axpy ( const Point &x, Factor &&... factor )
      {
        basisFunctionSet().axpy( x, std::forward< Factor >( factor )..., *this );
      }

      const BasisFunctionSetType &basisFunctionSet () const { return localMatrix_.rangeBasisFunctionSet(); }

    private:
      LocalMatrixType &localMatrix_;
      unsigned int col_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_COMMON_LOCALMATRIXCOLUMN_HH
