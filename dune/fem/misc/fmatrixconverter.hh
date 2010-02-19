#ifndef DUNE_FIELDMATRIXCONVERTER_HH
#define DUNE_FIELDMATRIXCONVERTER_HH

#include <cassert>
#include <dune/common/fmatrix.hh>

namespace Dune {

template <class VectorType, class ConvertToType> 
class FieldMatrixConverter;

//! convert a FieldVector with length n * m to a FieldMatrix with n rows and m cols 
template <typename K, int n, int m> 
class FieldMatrixConverter<FieldVector<K ,n * m> , FieldMatrix<K ,n, m> >
{
public:
  typedef FieldVector<K ,n * m> InteralVectorType;

  typedef InteralVectorType row_type;
  typedef row_type RowType;

  //! export the type representing the field
  typedef K field_type;

  //! export the type representing the components
  typedef K block_type;

  //! The type used for the index access and size operations.
  typedef std::size_t size_type;

  enum {
    //! The number of rows.
    rows = n,
    //! The number of columns.
    cols = m,
    //! The total dimension 
    dimension = n * m 
  };

  FieldMatrixConverter(InteralVectorType& v) 
    : vec_( v )
#ifndef NDEBUG 
    , mutableVec_( true )
#endif
  {
  }

  FieldMatrixConverter(const InteralVectorType& v) 
    : vec_(const_cast<InteralVectorType&> (v) )
#ifndef NDEBUG 
    , mutableVec_( false )
#endif
  {
  }

  FieldMatrixConverter(const FieldMatrixConverter& other) 
    : vec_( other.vec_ )
#ifndef NDEBUG 
    , mutableVec_( other.mutableVec_ )
#endif
  {}

  // return row 
  K* operator [] (const size_t row) 
  {
    assert( mutableVec_ );
    return &vec_[ row * cols ] ;
  }
  
  // return row  
  const K* operator [] (const size_t row) const 
  {
    return &vec_[ row * cols ];
  }

  FieldMatrixConverter& operator = (const FieldMatrixConverter& other) 
  {
    assert( mutableVec_ );
    vec_ = other.vec_;
#ifndef NDEBUG 
    mutableVec_ = other.mutableVec_;
#endif
  }

  FieldMatrixConverter& operator = (const FieldMatrix<K, n, m>& matrix) 
  {
    assert( mutableVec_ );
    for(size_t i=0; i< rows; ++i) 
    {
      for(size_t j=0; j<cols; ++j)
      {
        vec_[ i * cols + j ] = matrix[ i ][ j ];
      }
    }
  }

  /** \brief Sends the matrix to an output stream */
  friend std::ostream& operator<< (std::ostream& s, 
                                   const FieldMatrixConverter& a)
  {
    for (size_type i=0; i<n; ++i)
    {
      for(size_type j=0; j<m; ++j)
        s << a[ i ][ j ] << " "; 
      s << std::endl;
    }
    return s;
  }

protected:
  InteralVectorType& vec_;
#ifndef NDEBUG 
  bool mutableVec_;
#endif
};

// method that implements the operator= of a FieldMatrix taking a FieldMatrixConverter
template<class K, int n, int m>
inline void istl_assign_to_fmatrix(FieldMatrix<K,n,m>& A, 
    const FieldMatrixConverter< FieldVector<K, n * m> , FieldMatrix<K,n,m> >& B)
{
  for(size_t i=0; i<n; ++i) 
  {
    for(size_t j=0; j<m; ++j)
    {
      A[ i ][ j ] = B[ i ][ j ];
    }
  }
}



} // end namespace Dune 
#endif // end DUNE_FIELDMATRIXCONVERTER_HH
