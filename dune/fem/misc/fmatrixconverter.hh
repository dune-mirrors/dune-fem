#ifndef DUNE_FEM_FIELDMATRIXCONVERTER_HH
#define DUNE_FEM_FIELDMATRIXCONVERTER_HH

#include <cassert>

#include <dune/common/fmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class VectorType, class ConvertToType >
    class FieldMatrixConverter;

    template< class K, int m >
    class FieldMatrixConverterRow;

    template< class K, int n, int m >
    class FlatFieldMatrix;
  }


  template< class K, int n, int m >
  struct DenseMatVecTraits< Fem::FieldMatrixConverter< FieldVector< K, n *m >, FieldMatrix< K, n, m > > >
  {
    typedef Fem::FieldMatrixConverter< FieldVector< K, n *m >, FieldMatrix< K, n, m > > derived_type;
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
    typedef K *container_type;
    typedef K value_type;
    typedef size_t size_type;
  };

  template< class K, int n, int m >
  struct DenseMatVecTraits< Fem::FlatFieldMatrix< K, n, m > >
    : public DenseMatVecTraits< Fem::FieldMatrixConverter< FieldVector< K, n *m >, FieldMatrix< K, n, m > > >
  {
  };



  namespace Fem
  {

    // derive from dense vector to inherit functionality
    template< class K, int m >
    class FieldMatrixConverterRow
      : public Dune::DenseVector< FieldMatrixConverterRow< K, m > >
    {
      typedef FieldMatrixConverterRow< K, m > This;
      typedef Dune::DenseVector< FieldMatrixConverterRow< K, m > > Base;

    public:
      typedef typename Base::size_type size_type;

      using Base::operator=;
      using Base::operator*;

      FieldMatrixConverterRow ( K *ptr )
        : ptr_( ptr )
      {
        assert( ptr_ );
      }

      FieldMatrixConverterRow ( const This &other )
        : ptr_( other.ptr_ )
      {}

      FieldMatrixConverterRow &operator= ( const This &other )
      {
        for( size_t i = 0; i < size(); ++i )
          ptr_[ i ] = other[ i ];
        return *this;
      }

      template< class Impl >
      FieldMatrixConverterRow &operator= ( const DenseVector< Impl > &other )
      {
        assert( other.size() == size() );
        for( size_t i = 0; i < size(); ++i )
          ptr_[ i ] = other[ i ];
        return *this;
      }

      template< class V >
      K operator* ( const DenseVector< V > &x ) const
      {
        assert( size() == x.size() );
        K result( 0 );
        for( size_type i = 0; i < size(); ++i )
          result += ptr_[ i ] * x[ i ];
        return result;
      }

      // make this thing a vector
      size_t size () const { return m; }
      K &operator[] ( size_t i ) { assert( i < size() ); return ptr_[ i ]; }
      const K &operator[] ( size_t i ) const { assert( i < size() ); return ptr_[ i ]; }

    private:
      K *ptr_;
    };



    //! convert a FieldVector with length n * m to a FieldMatrix with n rows and m cols
    template< typename K, int n, int m >
    class FieldMatrixConverter< FieldVector< K, n *m >, FieldMatrix< K, n, m > >
      : public Dune::DenseMatrix< FieldMatrixConverter< FieldVector< K, n *m >, FieldMatrix< K, n, m > > >
    {
      typedef Dune::DenseMatrix< FieldMatrixConverter< FieldVector< K, n *m >, FieldMatrix< K, n, m > > > Base;

    public:
      //! internal storage of matrix
      typedef FieldVector< K, n *m > InteralVectorType;

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

      FieldMatrixConverter &operator= ( const FieldMatrixConverter &other )
      {
        assert( mutableVec_ );
        vec_ = other.vec_;
    #ifndef NDEBUG
        mutableVec_ = other.mutableVec_;
    #endif
        return *this;
      }

      FieldMatrixConverter &operator= ( const FieldMatrix< K, n, m > &matrix )
      {
        assert( mutableVec_ );
        for( int i = 0; i<rows; ++i )
        {
          for( int j = 0, ij = i * cols; j<cols; ++j, ++ij )
          {
            (*vec_)[ ij ] = matrix[ i ][ j ];
          }
        }

        return *this;
      }

      FieldMatrixConverter &operator+= ( const FieldMatrix< K, n, m > &matrix )
      {
        assert( mutableVec_ );
        for( int i = 0; i<rows; ++i )
        {
          for( int j = 0, ij = i * cols; j<cols; ++j, ++ij )
          {
            (*vec_)[ ij ] += matrix[ i ][ j ];
          }
        }
        return *this;
      }

      /** \brief Sends the matrix to an output stream */
      friend std::ostream &operator<< ( std::ostream &s,
                                        const FieldMatrixConverter &a )
      {
        for( size_type i = 0; i < n; ++i )
        {
          for( size_type j = 0; j < m; ++j )
            s << a[ i ][ j ] << " ";
          s << std::endl;
        }
        return s;
      }

      // make this thing a matrix
      size_type mat_rows () const { return rows; }
      size_type mat_cols () const { return cols; }

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

    template< typename K, int n, int m >
    class FlatFieldMatrix
      : public Dune::DenseMatrix< FlatFieldMatrix< K, n, m > >
    {
      typedef Dune::DenseMatrix< FlatFieldMatrix< K, n, m > > Base;

    public:
      //! internal storage of matrix
      typedef FieldVector< K, n * m > InteralVectorType;

      //! type of class return upon operator [] which behaves like a reference
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

    public:
      FlatFieldMatrix() {}
      FlatFieldMatrix( const K& val )
      {
        // all entries == val
        data_ = val;
      }

      FlatFieldMatrix( const FlatFieldMatrix& other )
      {
        data_ = other.data_;
      }

      FlatFieldMatrix( const FieldMatrix< K, n, m >& other )
      {
        assign( other );
      }

      FlatFieldMatrix &operator= ( const FlatFieldMatrix &other )
      {
        data_ = other.data_;
        return *this;
      }

      FlatFieldMatrix &operator= ( const K& val )
      {
        data_ = val;
        return *this;
      }

      FlatFieldMatrix &operator= ( const FieldMatrix< K, n, m > &other )
      {
        assign( other );
        return *this;
      }

      FlatFieldMatrix &operator+= ( const FieldMatrix< K, n, m > &other )
      {
        for( int i = 0; i<rows; ++i )
        {
          for( int j = 0, ij = i * cols; j<cols; ++j, ++ij )
            data_[ ij ] += other[ i ][ j ];
        }
        return *this;
      }

      void assign( const FieldMatrix< K, n, m > &other )
      {
        for( int i = 0; i<rows; ++i )
        {
          for( int j = 0, ij = i * cols; j<cols; ++j, ++ij )
            data_[ ij ] = other[ i ][ j ];
        }
      }

      K* data() { return &data_[ 0 ]; }
      const K* data() const { return &data_[ 0 ]; }

      // make this thing a matrix
      size_type mat_rows () const { return rows; }
      size_type mat_cols () const { return cols; }

      row_reference mat_access ( size_type i )
      {
        assert( i < rows );
        return row_reference( &(data_[ i * cols ]) );
      }

      const_row_reference mat_access ( size_type i ) const
      {
        assert( i < rows );
        return const_row_reference( &(data_[ i * cols ]) );
      }

      /** \brief Sends the matrix to an output stream */
      friend std::ostream &operator<< ( std::ostream &s,
                                        const FlatFieldMatrix &a )
      {
        for( size_type i = 0; i < n; ++i )
        {
          for( size_type j = 0; j < m; ++j )
            s << a[ i ][ j ] << " ";
          s << std::endl;
        }
        return s;
      }

    protected:
      InteralVectorType  data_;
    };


  } // namespace Fem


  namespace detail {

    template < class K, int n, int m, class Converter>
    inline void toMatrix( FieldMatrix< K, n, m>& A, const Converter& B )
    {
      for( int i = 0; i < n; ++i )
        for( int j = 0; j < m; ++j )
          A[ i ][ j ] = B[ i ][ j ];
    }
  }

  // method that implements the operator= of a FieldMatrix taking a FieldMatrixConverter
  template< class K, int n, int m >
  struct DenseMatrixAssigner< FieldMatrix< K, n, m >,
    Fem::FieldMatrixConverter< FieldVector< K, n*m >, FieldMatrix< K, n, m > > >
  {
    static void apply ( FieldMatrix< K, n, m > &A,
                        const Fem::FieldMatrixConverter< FieldVector< K, n *m >, FieldMatrix< K, n, m > > &B )
    {
      detail::toMatrix( A, B );
    }
  };

  template< class K, int n, int m >
  struct DenseMatrixAssigner< DenseMatrix< FieldMatrix< K, n, m > >,
    Fem::FieldMatrixConverter< FieldVector< K, n*m >, FieldMatrix< K, n, m > > >
  {
    static void apply ( DenseMatrix< FieldMatrix< K, n, m > > &A,
                        const Fem::FieldMatrixConverter< FieldVector< K, n *m >, FieldMatrix< K, n, m > > &B )
    {
      detail::toMatrix( A, B );
    }
  };

  // method that implements the operator= of a FieldMatrix taking a FlatFieldMatrix
  template< class K, int n, int m >
  struct DenseMatrixAssigner< FieldMatrix< K, n, m >, Fem::FlatFieldMatrix< K, n, m > >
  {
    static void apply ( FieldMatrix< K, n, m > &A,
                        const Fem::FlatFieldMatrix< K, n, m > &B )
    {
      detail::toMatrix( A, B );
    }
  };

  template< class K, int n, int m >
  struct DenseMatrixAssigner< DenseMatrix< FieldMatrix< K, n, m > >, Fem::FlatFieldMatrix< K, n, m > >
  {
    static void apply ( DenseMatrix< FieldMatrix< K, n, m > > &A,
                        const Fem::FlatFieldMatrix< K, n, m > &B )
    {
      detail::toMatrix( A, B );
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_FEM_FIELDMATRIXCONVERTER_HH
