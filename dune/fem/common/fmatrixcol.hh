#ifndef DUNE_FEM_COMMON_FMATRIXCOL_HH
#define DUNE_FEM_COMMON_FMATRIXCOL_HH

#include <type_traits>

#include <dune/common/densevector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class FieldMatrix >
  class FieldMatrixColumn;



  // DenseMatVecTraits for FieldMatrixColumn
  // ---------------------------------------

  template< class FieldMatrix >
  struct DenseMatVecTraits< FieldMatrixColumn< FieldMatrix > >
  {
    typedef FieldMatrixColumn< FieldMatrix > derived_type;

    typedef typename DenseMatVecTraits< typename std::remove_const< FieldMatrix >::type >::value_type value_type;
    typedef typename DenseMatVecTraits< typename std::remove_const< FieldMatrix >::type >::size_type size_type;
  };


  // FieldTraits for FieldMatrixColumn
  // ---------------------------------

  template< class FieldMatrix >
  struct FieldTraits< FieldMatrixColumn< FieldMatrix > >
  {
    typedef typename FieldTraits< typename std::remove_const< FieldMatrix >::type >::field_type field_type;
    typedef typename FieldTraits< typename std::remove_const< FieldMatrix >::type >::real_type real_type;
  };



  // FieldMatrixColumn
  // -----------------

  template< class K, int m, int n >
  class FieldMatrixColumn< FieldMatrix< K, m, n > >
    : public DenseVector< FieldMatrixColumn< FieldMatrix< K, m, n > > >
  {
    typedef DenseVector< FieldMatrixColumn< FieldMatrix< K, m, n > > > Base;

  public:
    static const int dimension = m;

    typedef typename Base::size_type size_type;
    typedef typename Base::value_type value_type;

    FieldMatrixColumn ( FieldMatrix< K, m, n > &fieldMatrix, int column )
      : fieldMatrix_( fieldMatrix ),
        column_( column )
    {}

    using Base::operator=;

    constexpr size_type size () const { return dimension; }

    const value_type &operator[] ( size_type i ) const { return fieldMatrix_[ i ][ column_ ]; }
    value_type &operator[] ( size_type i ) { return fieldMatrix_[ i ][ column_ ]; }

  private:
    FieldMatrix< K, m, n > &fieldMatrix_;
    int column_;
  };

  template< class K, int m, int n >
  class FieldMatrixColumn< const FieldMatrix< K, m, n > >
    : public DenseVector< FieldMatrixColumn< const FieldMatrix< K, m, n > > >
  {
    typedef DenseVector< FieldMatrixColumn< const FieldMatrix< K, m, n > > > Base;

  public:
    static const int dimension = m;

    typedef typename Base::size_type size_type;
    typedef typename Base::value_type value_type;

    FieldMatrixColumn ( const FieldMatrix< K, m, n > &fieldMatrix, int column )
      : fieldMatrix_( fieldMatrix ),
        column_( column )
    {}

    using Base::operator=;

    constexpr size_type size () const { return dimension; }
    const value_type &operator[] ( size_type i ) const { return fieldMatrix_[ i ][ column_ ]; }

  private:
    const FieldMatrix< K, m, n > &fieldMatrix_;
    int column_;
  };

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_FMATRIXCOL_HH
