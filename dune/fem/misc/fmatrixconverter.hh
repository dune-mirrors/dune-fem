#ifndef DUNE_FIELDMATRIXCONVERTER_HH
#define DUNE_FIELDMATRIXCONVERTER_HH

#include <dune/common/fmatrix.hh>

namespace Dune {

template <class VectorType, class ConvertToType> 
class FMatrixConverter;

//! convert a FieldVector with length n * m to a FieldMatrix with n rows and m cols 
template <typename K, int n, int m> 
class FMatrixConverter<FieldVector<K ,n * m> , FieldMatrix<K ,n, m> >
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

  FMatrixConverter(InteralVectorType& v) 
    : vec_( v )
#ifndef NDEBUG 
    , mutableVec_( true )
#endif
  {
  }

  FMatrixConverter(const InteralVectorType& v) 
    : vec_(const_cast<InteralVectorType&> (v) )
#ifndef NDEBUG 
    , mutableVec_( false )
#endif
  {
  }

  FMatrixConverter(const FMatrixConverter& org) 
    : vec_( org.vec_ )
#ifndef NDEBUG 
    , mutableVec_( org.mutableVec_ )
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

  FMatrixConverter& operator = (const FMatrixConverter& other) 
  {
    assert( mutableVec_ );
    vec_ = other.vec_;
#ifndef NDEBUG 
    mutableVec_ = org.mutableVec_;
#endif
  }

  FMatrixConverter& operator = (const FieldMatrix<K, n, m>& matrix) 
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

protected:
  InteralVectorType& vec_;
#ifndef NDEBUG 
  const bool constVec_;
#endif
};

// method that implements the operator= of a FieldMatrix taking a FMatrixConverter
template<class K, int n, int m>
inline void istl_assign_to_fmatrix(FieldMatrix<K,n,m>& A, 
    const FMatrixConverter< FieldVector<K, n * m> , FieldMatrix<K,n,m> >& B)
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
