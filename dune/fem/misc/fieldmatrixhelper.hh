#ifndef DUNE_FEM_FIELDMATRIXHELPER_HH
#define DUNE_FEM_FIELDMATRIXHELPER_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/dynmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    namespace FieldMatrixHelper
    {

      template< class Field1, class Field2, class Field3, int m, int n >
      inline void multiply ( const FieldMatrix< Field1, m, n > &A,
                             const FieldVector< Field2, n > &x,
                             FieldVector< Field3, m > &y )
      {
        for( int i = 0; i < m; ++i )
        {
          Field3 &value = y[ i ];

          value = 0;
          for( int j = 0; j < n; ++j )
            value += A[ i ][ j ] * x[ j ];
        }
      }



      template< class Field1, class Field2, class Field3, int m, int n, int p >
      inline void multiply ( const FieldMatrix< Field1, m, n > &A,
                             const FieldMatrix< Field2, n, p > &B,
                             FieldMatrix< Field3, m, p > &C )
      {
        for( int i = 0; i < m; ++i )
        {
          for( int j = 0; j < p; ++j )
          {
            Field3 &value = C[ i ][ j ];

            value = 0;
            for( int k = 0; k < n; ++k )
              value += A[ i ][ k ] * B[ k ][ j ];
          }
        }
      }

      template< class Field1, class Field2, class Field3 >
      inline void multiply ( const DynamicMatrix< Field1 >& A,
                             const DynamicMatrix< Field2 >& B,
                             DynamicMatrix< Field3 >& C )
      {
        const int m = A.rows();
        const int n = A.cols();
        const int p = B.cols();

        C.resize( m, p );

        assert( A.cols() == B.rows() );
        assert( A.rows() == C.rows() );
        assert( B.cols() == C.cols() );

        for( int i = 0; i < m; ++i )
        {
          for( int j = 0; j < p; ++j )
          {
            Field3 &value = C[ i ][ j ];

            value = 0;
            for( int k = 0; k < n; ++k )
              value += A[ i ][ k ] * B[ k ][ j ];
          }
        }
      }

      template< class Field1, class Field2, class Field3, int m, int n, int p >
      inline void multiply ( const FieldMatrix< Field1, m, n > &A,
                             const FieldMatrix< Field2, n, p > &B,
                             Field3* C)
      {
        for( int i = 0, ip = 0; i < m; ++i )
        {
          for( int j = 0; j < p; ++j , ++ ip )
          {
            Field3 &value = C[ ip ];
            value = 0;
            for( int k = 0; k < n; ++k )
              value += A[ i ][ k ] * B[ k ][ j ];
          }
        }
      }

    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FIELDMATRIXHELPER_HH
