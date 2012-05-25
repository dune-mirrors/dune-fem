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

template < class K > 
class SSEMatrix 
{
  SSEMatrix(const SSEMatrix& other) ;
public:
  //typedef K  Field __attribute__((aligned (16)));
  typedef K  Field ;
  typedef Field* RowType; 
  typedef RowType* MatrixType;

  typedef MatrixType SSEType;
protected:  
  typedef void multiplySSE_t(const size_t I, const size_t J, 
                             K** const &A, 
                             K*  const &X, 
                             K*  &Y);

  template <int n, int m> 
  struct Multiply
  {
    static multiplySSE_t* mv(const int rows, const int cols) 
    {
      if( n == rows && m == cols )
      {
        return multSSE;
      }
      else if( n == rows ) 
      {
        return Multiply< n, m - 1>::mv( rows, cols );
      }
      else 
      {
        return Multiply<n-1, m> ::mv( rows, cols );
      }
    }

    static multiplySSE_t* umtv(const int rows, const int cols) 
    {
      if( n == rows && m == cols )
      {
        return multUmtvSSE;
      }
      else if( n == rows ) 
      {
        return Multiply< n, m - 1>::umtv( rows, cols );
      }
      else 
      {
        return Multiply<n-1, m> ::umtv( rows, cols );
      }
    }

    static inline 
    void multSSE(const size_t I, const size_t J, 
                 K** const &A, 
                 K*  const &X, 
                 K*  &Y)
    {
      assert( n == I );
      assert( m == J );
      SSEMatrix< K > :: multiplySSE(n, m, A, X , Y);
    }

    static inline 
    void multUmtvSSE(const size_t I, const size_t J, 
                     K** const &A, 
                     K*  const &X, 
                     K*  &Y)
    {
      assert( n == I );
      assert( m == J );
      SSEMatrix< K > :: umtv(n, m, A, X , Y);
    }
  };

  template <int m> 
  struct Multiply<0, m> 
  {
    static multiplySSE_t* mv(const int rows, const int cols) 
    {
      return 0; 
    }
    static multiplySSE_t* umtv(const int rows, const int cols) 
    {
      return 0; 
    }
  };

  template <int n> 
  struct Multiply<n, 0> 
  {
    static multiplySSE_t* mv(const int rows, const int cols) 
    {
      return 0; 
    }
    static multiplySSE_t* umtv(const int rows, const int cols) 
    {
      return 0; 
    }
  };

  MatrixType sseMat_;
  mutable SSEVector< K > xTmp_;
  mutable SSEVector< K > yTmp_;

  mutable MatrixType quadMat_;

  multiplySSE_t* mvSSE_;
  multiplySSE_t* umtvSSE_;

  const size_t rows_;
  const size_t cols_;
  const size_t realRows_;

  enum { countSize = 100 };
public:
  SSEMatrix(const size_t n, const size_t m ,
            const size_t realRows ) 
    : xTmp_(  m ) 
    , yTmp_(  n )
    , quadMat_( 0 )
    , mvSSE_(   Multiply<countSize, countSize>::mv( realRows, m ) )
    , umtvSSE_( Multiply<countSize, countSize>::umtv( realRows, m ) )
    , rows_( n ), cols_ ( m ), realRows_( realRows )
  {
    assert( mvSSE_ );
    sseMat_  = new RowType [ n ];
    quadMat_ = new RowType [ n ];
    assert( sseMat_ );
    const size_t SIZE_J = (2*(( m + 1 ) / 2));
    for( size_t i = 0; i < rows_; ++i) 
    {
      sseMat_[ i ] = new K [ SIZE_J ];
      assert( sseMat_[ i ] );
    }
  } 

  // dummy method for fillRangeCache
  void resize( const size_t ) {}

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

  template <class X, class Y> 
  void mv(const X& x, Y& y, const size_t offset ) const 
  {
    mv( sseMat_, rows(), x, y, offset );
  }

  template <class X, class Y, class quad_t > 
  void mv(const X& x, Y& y, 
          const quad_t& quad,
          const size_t offset ) const 
  {
    assert( realRows_ == quad.nop() );
    for( size_t i = 0; i < realRows_; ++i)
    {
      quadMat_[ i ] = sseMat_[ quad.cachingPoint( i ) ];
    }
    mv( quadMat_, realRows_, x, y, offset );
  }

#if 1
  template <class X, class Y> 
  void 
  mv(const MatrixType& mat, 
     const size_t realRows,
     const X& x, Y& y, 
     const size_t offset ) const 
  {
    for(size_t r = 0; r < offset ; ++r )
    {
      for(size_t c = 0, cR = 0; c < cols(); ++c, cR += r ) xTmp_[ c ] = x[ cR ];
      //multiplySSE( realRows, cols(), mat, xTmp_.raw(), yTmp_.raw() );
      mvSSE_( realRows, cols(), mat, xTmp_.raw(), yTmp_.raw() );
      for(size_t i = 0; i < realRows; ++i ) y[ i ][ r ] = yTmp_[ i ];
    } 
  }
#else
  template <int n, int m, class X, class Y> 
  void 
  mv(const MatrixType& mat, 
     const size_t realRows,
     const X& x, Y& y, 
     const size_t offset ) const 
  {
    for (size_t i = 0; i < realRows; ++i)
    {
      const RowType& matRow = mat[ i ];

      for(size_t r = 0; r < offset; ++r ) 
      {
        y[ i ][ r ] = 0;
      }

      for (size_t j = 0; j < cols_; ++j)
      {
        for(size_t r = 0; r < offset; ++r ) 
        {
          y[ i ][ r ] += matRow[ j ] * x[ j * offset + r ];
        }
      }
    }
  }
#endif

  template <class X, class Y, class quad_t > 
  void umtv(const X& x, Y& y, const quad_t& quad, 
            const size_t offset ) const 
  {
    const size_t quadSize = quad.nop();
    for( size_t i = 0; i < quadSize; ++i) 
    {
      quadMat_[ i ] = sseMat_[ quad.cachingPoint( i ) ];
    }
    umtv( quadMat_, quadSize, x, y, offset );
  }

  template <class X, class Y, class quad_t > 
  void umtv(const X& x, Y& y, const quad_t& quad, 
            const size_t offset , bool ) const 
  {
    const size_t quadSize = quad.nop();
    for( size_t i = 0; i < quadSize; ++i) 
    {
      quadMat_[ i ] = sseMat_[ quad.cachingPoint( i ) ];
    }
    umtv( quadMat_, quadSize, x, y, offset , true );
  }

  template <class X, class Y> 
  void umtv(const X& x, Y& y, const size_t offset ) const 
  {
    umtv( sseMat_, rows(), x, y, offset );
  }

  template <class X, class Y> 
  void umtv(const X& x, Y& y, const size_t offset , bool ) const 
  {
    umtv( sseMat_, rows(), x, y, offset , true  );
  }

#if 1
  template <class X, class Y> 
  void umtv(const MatrixType& mat, const size_t realRows,
            const X& x, Y& y, const size_t offset , bool ) const 
  {
    for(size_t r = 0; r < offset ; ++r )
    {
      for(size_t i = 0, iR = 0; i < realRows; ++i, iR += r ) yTmp_[ i ] = x[ iR ];
      //multiplySSE( realRows, cols(), mat, xTmp_.raw(), yTmp_.raw() );
      umtvSSE_( realRows, cols(), mat, yTmp_.raw(), xTmp_.raw() );
      for(size_t c = 0, cR = 0; c < cols(); ++c, cR += r ) y[ cR ] += xTmp_[ c ];
    } 
  }

  template <class X, class Y> 
  void umtv(const MatrixType& mat, const size_t realRows,
            const X& x, Y& y, const size_t offset ) const 
  {
    for(size_t r = 0; r < offset ; ++r )
    {
      for(size_t i = 0; i < realRows; ++i ) yTmp_[ i ] = x[ i ][ r ];
      //multiplySSE( realRows, cols(), mat, xTmp_.raw(), yTmp_.raw() );
      umtvSSE_( realRows, cols(), mat, yTmp_.raw(), xTmp_.raw() );
      for(size_t c = 0, cR = 0; c < cols(); ++c, cR += r ) y[ cR ] += xTmp_[ c ];
    } 
  }

  static inline 
  void umtv(const size_t rows,
            const size_t cols,
            K** const &A, 
            K*  const &X, 
            K*  &Y)
  {
    for( size_t j = 0; j < cols; ++j ) 
    {
      Y[ j ] = 0;
    }
    for( size_t i = 0; i < rows; ++i )
    {
      for( size_t j = 0; j < cols; ++j ) 
      {
        Y[ j ] += A[ i ][ j ] * X[ i ];
      }
    }
  }
#else 
  template <class X, class Y> 
  void umtv(K** const & mat, const size_t realRows,
            const X& x, Y& y, const size_t offset ) const 
  {
    for( size_t row = 0; row < realRows ; ++row )
    {
      for( size_t col = 0; col < cols_; ++col ) 
      {
        const Field& matvalue = mat[ row ][ col ];
        size_t colR = col * offset;
        for( int r = 0; r < offset ; ++r, ++ colR ) 
        {
          y[ colR ] += matvalue * x[ row ][ r ];
        }
      }
    }
  }
#endif

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

  static inline 
  void multiplySSE(const size_t I, const size_t J, 
                   double** const &A, 
                   double*  const &X, 
                   double*  &Y)
  {
#if 1
    for(size_t i = 0; i < I; ++i ) 
    {
      Y[ i ] = 0; 
      for( size_t j = 0; j < J; ++j ) 
      {
        Y[ i ] += A[ i ][ j ] * X[ j ]; 
      }
    }
#else 
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
#endif
  }

  static inline 
  void multiplySSE(const size_t I, const size_t J, 
                   float** const &A, 
                   float* const &X, 
                   float* Y)
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

};

} // end namespace Fem 

} // end namespace Dune 
#endif
