#ifndef DUNE_FEM_OPERATOR_COMMON_LOCALMATRIXCOLUMN_HH
#define DUNE_FEM_OPERATOR_COMMON_LOCALMATRIXCOLUMN_HH

#include <utility>

#include <dune/fem/operator/common/temporarylocalmatrix.hh>

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

      typedef RangeFieldType value_type;
      typedef unsigned int size_type;

      LocalMatrixColumn ( LocalMatrixType &localMatrix, unsigned int col )
        : localMatrix_( localMatrix ), col_( col )
      {}

      RangeFieldType operator[] ( size_type row ) const { return localMatrix_.get( row, col_ ); }

      LocalMatrixEntry< LocalMatrixType > operator[] ( size_type row ) { return LocalMatrixEntry< LocalMatrixType >( localMatrix_, row, col_ ); }

      template< class Point, class... Factor >
      void axpy ( const Point &x, Factor &&... factor )
      {
        basisFunctionSet().axpy( x, std::forward< Factor >( factor )..., *this );
      }

      template< class Quadrature, class... Factor >
      void axpyQuadrature ( const Quadrature &quadrature, Factor &&... factor )
      {
        basisFunctionSet().axpy( quadrature, std::forward< Factor >( factor )..., *this );
      }

      const BasisFunctionSetType &basisFunctionSet () const { return localMatrix_.rangeBasisFunctionSet(); }

      auto size() const { return basisFunctionSet().size(); }

    private:
      LocalMatrixType &localMatrix_;
      unsigned int col_;
    };



    // LocalMatrixColumn for TemporaryLocalMatrix
    // ------------------------------------------

    template< class DomainSpace, class RangeSpace >
    class LocalMatrixColumn< TemporaryLocalMatrix< DomainSpace, RangeSpace > >
    {
      typedef LocalMatrixColumn< TemporaryLocalMatrix< DomainSpace, RangeSpace > > ThisType;

    public:
      typedef TemporaryLocalMatrix< DomainSpace, RangeSpace > LocalMatrixType;

      typedef typename LocalMatrixType::RangeFieldType RangeFieldType;
      typedef typename LocalMatrixType::RangeBasisFunctionSetType BasisFunctionSetType;

      typedef typename BasisFunctionSetType::RangeType RangeType;
      typedef typename BasisFunctionSetType::JacobianRangeType JacobianRangeType;
      typedef typename BasisFunctionSetType::HessianRangeType HessianRangeType;

      typedef RangeFieldType value_type;
      typedef unsigned int size_type;

      LocalMatrixColumn ( LocalMatrixType &localMatrix, unsigned int col )
        : localMatrix_( localMatrix ), col_( col )
      {}

      const value_type &operator[] ( size_type row ) const { return localMatrix_[ row ][ col_ ]; }
      value_type &operator[] ( size_type row ) { return localMatrix_[ row ][ col_ ]; }

      template< class Point, class... Factor >
      void axpy ( const Point &x, Factor &&... factor )
      {
        basisFunctionSet().axpy( x, std::forward< Factor >( factor )..., *this );
      }

      template< class Quadrature, class... Factor >
      void axpyQuadrature ( const Quadrature &quadrature, Factor &&... factor )
      {
        basisFunctionSet().axpy( quadrature, std::forward< Factor >( factor )..., *this );
      }

      const BasisFunctionSetType &basisFunctionSet () const { return localMatrix_.rangeBasisFunctionSet(); }

      auto size() const { return basisFunctionSet().size(); }

    private:
      LocalMatrixType &localMatrix_;
      unsigned int col_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_COMMON_LOCALMATRIXCOLUMN_HH
