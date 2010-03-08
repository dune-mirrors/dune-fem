#ifndef DUNEFEM_SSEMATRIXVECTOR_HH
#define DUNEFEM_SSEMATRIXVECTOR_HH

#include <iostream>
#include <cassert>

extern "C" {
  #include <emmintrin.h>
  #include <mmintrin.h>
}

namespace Dune { 

namespace Fem {  

template <class K> 
class SSEVector
{
  SSEVector(const SSEVector& );
public:
  //typedef K Field __attribute__((aligned (16)));
  typedef K Field;
  typedef Field* RowType;

  typedef RowType SSEType;
protected:  
  RowType sseVec_;
  const size_t size_;
public:
  SSEVector(const size_t SIZE) : size_( SIZE ) 
  {
    const size_t SIZE_J = (2*((SIZE+1)/2));
    sseVec_ = new K [ SIZE_J ];
    assert( sseVec_ );
  } 

  ~SSEVector () 
  { 
    delete [] sseVec_;
  } 

  void resize( const size_t ) {}

  K& operator [] (const size_t i) 
  {
    assert( i < size() );
    return sseVec_[ i ];
  }

  const K& operator [] (const size_t i) const 
  {
    assert( i < size() );
    return sseVec_[ i ];
  }

  RowType& raw() { return sseVec_; }
  const RowType& raw() const { return sseVec_; }
  size_t size() const { return size_; }
};

template<class RowVector, class ColVector>  
inline 
void multiplySSE(const size_t I, const size_t J, 
                 float** const &A, 
                 const RowVector& X, 
                 ColVector& Y)
{
  const size_t N = 3;

  const size_t FI = ((N)*((I)/(N)));
  const size_t FJ = (4*((J)/4));

  const size_t AI = ((N)*((I+(N-1)/(N))));
  const size_t AJ = (4*((J+3)/4));

  for (size_t i = 0; i < FI; i += N) 
  {
    __m128 tmp[N];
    for (size_t n = 0; n < N; n++)
      tmp[n] = _mm_set1_ps(0);
      
    for (size_t j = 0; j < FJ; j += 4)
      for (size_t n = 0; n < N; n++)
        tmp[n] = _mm_add_ps(tmp[n], _mm_mul_ps(_mm_load_ps(A[i+n]+j), _mm_load_ps(X+j)));
        
    for (size_t n = 0; n < N; n++)
      _mm_store_ss(Y+i+n, _mm_add_ss(_mm_add_ss(_mm_add_ss(tmp[n], _mm_shuffle_ps(tmp[n], tmp[n], 1)), 
                  _mm_shuffle_ps(tmp[n], tmp[n], 2)), _mm_shuffle_ps(tmp[n], tmp[n], 3)));

    for (size_t n = 0; n < N; n++)
      for (size_t j = FJ; j < J; j++)
        Y[i+n] += A[i+n][j] * X[j];
  }
  for (size_t i = FI; i < I; i++) {
    __m128 tmp = _mm_set1_ps(0);

    for (size_t j = 0; j < FJ; j += 4)
      tmp = _mm_add_ps(tmp, _mm_mul_ps(_mm_load_ps(A[i]+j), _mm_load_ps(X+j)));

    _mm_store_ss(Y+i, _mm_add_ss(_mm_add_ss(_mm_add_ss(tmp, _mm_shuffle_ps(tmp, tmp, 1)), _mm_shuffle_ps(tmp, tmp, 2)), _mm_shuffle_ps(tmp, tmp, 3)));

    for (size_t j = FJ; j < J; j++)
      Y[i] += A[i][j] * X[j];
  }
}

template<class RowVector, class ColVector>  
inline 
void multiplySSE(const size_t I, const size_t J, 
                 double** const &A, 
                 const RowVector&X, 
                 ColVector& Y)
{
  const size_t N = 3;
  const size_t FI = ((N)*(( I )/(N)));  
  const size_t FJ = (2*( J )/2);

  for (size_t i = 0; i < FI; i += N) 
  {
    __m128d tmp[N];
    for (size_t n = 0; n < N; ++n)
      tmp[n] = _mm_set1_pd(0);

    for (size_t j = 0; j < FJ; j += 2)
      for (size_t n = 0; n < N; ++n)
        tmp[n] = _mm_add_pd(tmp[n], _mm_mul_pd(_mm_load_pd(A[i+n]+j), _mm_load_pd(X+j)));

    for (size_t n = 0; n < N; ++n)
      _mm_store_sd(Y+i+n, _mm_add_sd(tmp[n], _mm_shuffle_pd(tmp[n], tmp[n], 1)));

    for (size_t n = 0; n < N; ++n)
      for (size_t j = FJ; j < J; ++j)
        Y[i+n] += A[i+n][j] * X[j];
  }

  for (size_t i = FI; i < I; ++i) 
  {
    __m128d tmp = _mm_set1_pd(0);

    for (size_t j = 0; j < FJ; j += 2)
      tmp = _mm_add_pd(tmp, _mm_mul_pd(_mm_load_pd(A[i]+j), _mm_load_pd(X+j)));

    _mm_store_sd(Y+i, _mm_add_sd(tmp, _mm_shuffle_pd(tmp, tmp, 1)));

    for (size_t j = FJ; j < J; ++j)
      Y[i] += A[i][j] * X[j];
  }
}

template<class Matrix, class RowVector, class ColVector>  
inline 
void multiplySSE(const Matrix& mat, 
                 const RowVector& X, 
                 ColVector& Y)
{
  multiplySSE( mat.rows(), mat.cols(),
               mat.raw(), X, Y);
}

template < class K > 
class SSEMatrix 
{
public:
  //typedef K  Field __attribute__((aligned (16)));
  typedef K  Field ;
  typedef Field* RowType; 
  typedef RowType* MatrixType;

  typedef MatrixType SSEType;
protected:  
  MatrixType sseMat_;
  mutable SSEVector< K >& xTmp_;
  mutable SSEVector< K >& yTmp_;

  mutable MatrixType quadMat_;

  const size_t rows_;
  const size_t cols_;

  mutable size_t qSize_;

public:

  SSEMatrix(const size_t n, const size_t m ) 
    : xTmp_(*(new SSEVector< K >( m ))) 
    , yTmp_(*(new SSEVector< K >( n )))
    , quadMat_( 0 ), rows_( n ), cols_ ( m ), qSize_( 0 )
  {
    sseMat_ = new RowType [ n ];
    assert( sseMat_ );
    const size_t SIZE_J = (2*(( m + 1 ) / 2));
    for( size_t i = 0; i < rows_; ++i) 
    {
      sseMat_[ i ] = new K [ SIZE_J ];
      assert( sseMat_[ i ] );
    }
  } 

  SSEMatrix(const size_t n, const size_t m,
            MatrixType sseMat, 
            SSEVector< K >& xTmp, 
            SSEVector< K >& yTmp)
    : sseMat_( sseMat ) 
    , xTmp_( xTmp ), yTmp_( yTmp )
    , rows_( n ), cols_ ( m ), quadMat_( 0 ), qSize_ ( 0 )
  {
  }

  SSEMatrix(const SSEMatrix& other) 
    : sseMat_( other.sseMat_ )
    , xTmp_( other.xTmp_ )
    , yTmp_( other.yTmp_ )
    , quadMat_( 0 )
    , rows_( other.rows_ )
    , cols_( other.cols_ )
    , qSize_( 0 )
  {}

  void resize( const size_t ) {}

  template <class quad_t>
  const SSEMatrix<K> getQuadMat( const quad_t& quad ) const 
  {
    if( ! quadMat_ ) 
    { 
      qSize_ = quad.nop(); 
      quadMat_ = new RowType [ qSize_ ];
    }

    assert( quadMat_ ); 
    for( size_t i = 0; i < qSize_; ++i) 
    {
      quadMat_[ i ]  = sseMat_[ quad.cachingPoint( i ) ];
    }
    return SSEMatrix<K>(qSize_, cols(), quadMat_, xTmp_, yTmp_);
  }

  ~SSEMatrix() 
  {
    if( quadMat_ ) 
    {
      for( size_t i = 0; i < rows_; ++i)
      {
        delete [] sseMat_[ i ];
      }
      delete [] sseMat_;
      delete [] quadMat_;

      delete & xTmp_; 
      delete & yTmp_;
    }
    sseMat_ = 0;
    quadMat_ = 0;
  }

  RowType& operator [] (const size_t i) 
  {
    assert( i < rows_ );
    return sseMat_[ i ];
  }

  const RowType& operator [] (const size_t i) const 
  {
    assert( i < cols_ );
    return sseMat_[ i ];
  }

  size_t rows () const { return rows_; }
  size_t cols () const { return cols_; }

#if 1
  template <class X, class Y> 
  void mv(const X& x, Y& y, const size_t offset ) const 
  {
    for(size_t r = 0; r < offset ; ++r )
    {
      for(size_t c = 0; c < cols(); ++c ) xTmp_[ c ] = x[ c * offset + r ];
      multiplySSE( rows(), cols(), raw(), xTmp_.raw(), yTmp_.raw() );
      for(size_t i = 0; i < rows(); ++i ) y[ i ][ r ] = yTmp_[ i ];
    } 
  }

#else
  template <class X, class Y> 
  void mv(const X& x, Y& y, const size_t offset ) const 
  {
    for (size_t i = 0; i < rows(); ++i)
    {
      const RowType& matRow = sseMat_[ i ];

      for(size_t r = 0; r < offset; ++r ) 
      {
        y[ i ][ r ] = 0;
      }

      for (size_t j = 0; j < cols(); ++j)
      {
        for(size_t r = 0; r < offset; ++r ) 
        {
          y[ i ][ r ] += matRow[ j ] * x[ j * offset + r ];
        }
      }
    }
  }
#endif

  template <class X, class Y> 
  void umtv(const X& x, Y& y, const size_t offset ) const 
  {
    for( size_t row = 0; row < rows() ; ++row )
    {
      for( size_t col = 0; col < cols(); ++col ) 
      {
        const Field& matvalue = sseMat_[ row ][ col ];
        for( int r = 0; r < offset ; ++r ) 
        {
          y[ col * offset + r ] += matvalue * x[ row ][ r ];
        }
      }
    }
  }

  template <class X, class Y> 
  void mv(const X& x, Y& y) const 
  {
    mv(x, y, 1);
  }

  template <class Y> 
  void mv(const SSEVector<K>& x, Y& y) const 
  {
    multiplySSE( rows(), cols(), raw(), x.raw(), y );
  }

  MatrixType& raw() { return sseMat_; }
  const MatrixType& raw() const { return sseMat_; }

  size_t size() const { return rows() * cols(); }
};

} // end namespace Fem 

} // end namespace Dune 
#endif
