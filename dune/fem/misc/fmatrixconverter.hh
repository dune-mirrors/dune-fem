#ifndef DUNE_FEM_FIELDMATRIXCONVERTER_HH
#define DUNE_FEM_FIELDMATRIXCONVERTER_HH

#include <cassert>

#include <dune/common/fmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    // forward declaration
    template< class VectorType, class ConvertToType >
    class FieldMatrixConverter;

    template< class K, int m >
    class FieldMatrixConverterRow;
  }
  
  
  template< class K, int n, int m >
  struct DenseMatVecTraits< Fem::FieldMatrixConverter< FieldVector< K, n*m >, FieldMatrix< K, n, m > > >
  {
    typedef Fem::FieldMatrixConverter< FieldVector< K, n*m >, FieldMatrix< K, n, m > > derived_type;
    typedef DenseMatVecTraits< FieldMatrix< K, n, m > > Traits;
    typedef typename Traits::container_type container_type;
    typedef typename Traits::value_type value_type;
    typedef typename Traits::size_type size_type;
    typedef typename Traits::row_type row_type;

    typedef Fem::FieldMatrixConverterRow< K, m > row_reference;
    typedef Fem::FieldMatrixConverterRow< K, m > const_row_reference;
  };

  template< class K, int m >
  struct DenseMatVecTraits< Fem::FieldMatrixConverterRow< K, m > >
  {
    typedef Fem::FieldMatrixConverterRow< K, m > derived_type;
    typedef K* container_type;
    typedef K value_type;
    typedef size_t size_type;
  };  


  namespace Fem
  {

    // derive from dense vector to inherit functionality 
    template< class K, int m >
    class FieldMatrixConverterRow
    : public DenseVector< FieldMatrixConverterRow< K , m > >
    {
      typedef DenseVector< FieldMatrixConverterRow< K , m > > Base;

    public:  
      typedef typename Base::size_type size_type;

      using Base::operator =;
      using Base::operator *;

      FieldMatrixConverterRow ( K *ptr )
      : ptr_( ptr )
      {
        assert( ptr_ );
      }

      FieldMatrixConverterRow& operator = ( const FieldMatrixConverterRow& other )
      {
        for(size_t i=0; i<vec_size(); ++i ) 
        {
          vec_access( i ) = other[ i ];
        }
        return *this;
      }

      template <class Impl> 
      FieldMatrixConverterRow& operator = ( const DenseVector< Impl >& other )
      {
        assert( other.size() == vec_size() );
        for(size_t i=0; i<vec_size(); ++i ) 
        {
          vec_access( i ) = other[ i ];
        }
        return *this;
      }

      template <class V>
      K operator* ( const DenseVector< V >& x ) const 
      {
        assert( vec_size() == x.size() );
        K result( 0 );
        for(size_type i=0; i<vec_size(); ++i )
        {
          result += vec_access( i ) * x[ i ];
        }
        return result;
      }

      // make this thing a vector
      size_t vec_size() const { return m; }
      K & vec_access(size_t i) { assert( i < vec_size() ); return ptr_[ i ]; }
      const K & vec_access(size_t i) const { assert( i < vec_size() ); return ptr_[ i ]; }

    private:
      K *ptr_;
    };




    //! convert a FieldVector with length n * m to a FieldMatrix with n rows and m cols 
    template< typename K, int n, int m > 
    class FieldMatrixConverter< FieldVector< K, n*m >, FieldMatrix< K, n, m > >
    : public DenseMatrix< FieldMatrixConverter< FieldVector< K , n*m >, FieldMatrix< K, n, m > > >
    {
      typedef DenseMatrix< FieldMatrixConverter< FieldVector< K , n*m >, FieldMatrix< K, n, m > > > Base;

    public:
      //! internal storage of matrix 
      typedef FieldVector<K ,n * m> InteralVectorType;

      //! type of class return upon operator [] which behaves like a reference
      //typedef FieldMatrixConverterRow< K , m> RowType;

      typedef typename Base::row_type row_type;
      typedef typename Base::row_reference row_reference;
      typedef typename Base::const_row_reference const_row_reference;

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

      FieldMatrixConverter ( InteralVectorType &v )
      : vec_( &v )
#ifndef NDEBUG 
        , mutableVec_( true )
#endif
      {}

      FieldMatrixConverter ( const InteralVectorType &v )
      : vec_( const_cast< InteralVectorType * >( &v ) )
#ifndef NDEBUG 
        , mutableVec_( false )
#endif
      {}

      FieldMatrixConverter ( const FieldMatrixConverter &other )
      : vec_( other.vec_ )
#ifndef NDEBUG 
        , mutableVec_( other.mutableVec_ )
#endif
      {}

#if 0
      // return row 
      RowType operator [] (const size_t row) 
      {
        assert( mutableVec_ );
        return RowType( (&vec_[ row * cols ]) );
      }
      
      // return row  
      RowType operator [] (const size_t row) const 
      {
        return RowType( (&vec_[ row * cols ]) );
      }

      template<class X, class Y>
      void mv(const X& x, Y& y) const 
      {
        for(size_type i=0; i<n; ++i) 
        {
          y[ i ] = (*this)[ i ] * x;
        }
      }

      template<class X, class Y>
      void umv(const X& x, Y& y) const 
      {
        for(size_type i=0; i<n; ++i) 
        {
          y[ i ] += (*this)[ i ] * x;
        }
      }

      template<class X, class Y>
      void mtv(const X& x, Y& y) const 
      {
        for( size_type i = 0; i < cols; ++i )
        {
          y[ i ] = 0;
          for( size_type j = 0; j < rows; ++j )
            y[ i ] += (*this)[ j ][ i ] * x[ j ];
        }
      }

      template<class X, class Y>
      void umtv(const X& x, Y& y) const 
      {
        for( size_type i = 0; i < cols; ++i )
        {
          for( size_type j = 0; j < rows; ++j )
            y[ i ] += (*this)[ j ][ i ] * x[ j ];
        }
      }

      K frobenius_norm() const
      {
        K ret(0);
        for( size_type i = 0; i < cols; ++i )
        {
          for( size_type j = 0; j < rows; ++j )
            ret += (*this)[ i ][ j ] * (*this)[ i ][ j ];
        }
        return std::sqrt(ret);
      }
#endif 

      FieldMatrixConverter& operator = (const FieldMatrixConverter& other) 
      {
        assert( mutableVec_ );
        vec_ = other.vec_;
    #ifndef NDEBUG 
        mutableVec_ = other.mutableVec_;
    #endif
        return *this;
      }

      FieldMatrixConverter& operator = (const FieldMatrix<K, n, m>& matrix) 
      {
        assert( mutableVec_ );
        for(size_t i=0; i< rows; ++i) 
        {
          for(size_t j=0; j<cols; ++j)
          {
            (*vec_)[ i * cols + j ] = matrix[ i ][ j ];
          }
        }
        return *this;
      }
      FieldMatrixConverter& operator += (const FieldMatrix<K, n, m>& matrix) 
      {
        assert( mutableVec_ );
        for(size_t i=0; i< rows; ++i) 
        {
          for(size_t j=0; j<cols; ++j)
          {
            (*vec_)[ i * cols + j ] += matrix[ i ][ j ];
          }
        }
        return *this;
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

      // make this thing a matrix
      size_type mat_rows() const { return rows; }
      size_type mat_cols() const { return cols; }

      row_reference mat_access ( size_type i )
      {
        assert( i < rows );
        return row_reference( (&(*vec_)[ i * cols ]) );
      }

      const_row_reference mat_access ( size_type i ) const
      {
        assert( i < rows );
        return const_row_reference( (&(*vec_)[ i * cols ]) );
      }

    protected:
      InteralVectorType *vec_;
    #ifndef NDEBUG 
      bool mutableVec_;
    #endif
    };

  } // namespace Fem

  // method that implements the operator= of a FieldMatrix taking a FieldMatrixConverter
  template<class K, int n, int m>
  inline void istl_assign_to_fmatrix(FieldMatrix<K,n,m>& A, 
      const Fem::FieldMatrixConverter< FieldVector<K, n * m> , FieldMatrix<K,n,m> >& B)
  {
    for(size_t i=0; i<n; ++i) 
    {
      for(size_t j=0; j<m; ++j)
      {
        A[ i ][ j ] = B[ i ][ j ];
      }
    }
  }
  template<class K, int n, int m>
  inline void istl_assign_to_fmatrix(DenseMatrix<FieldMatrix<K,n,m> >& A, 
      const Fem::FieldMatrixConverter< FieldVector<K, n * m> , FieldMatrix<K,n,m> >& B)
  {
    for(size_t i=0; i<n; ++i) 
    {
      for(size_t j=0; j<m; ++j)
      {
        A[ i ][ j ] = B[ i ][ j ];
      }
    }
  }

} // namespace Dune

#endif // #ifndef DUNE_FEM_FIELDMATRIXCONVERTER_HH
