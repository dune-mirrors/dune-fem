#ifndef DUNE_FEM_FIELDMATRIXHELPER_HH
#define DUNE_FEM_FIELDMATRIXHELPER_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune
{

  namespace FieldMatrixHelper
  {

    template< class Field, int m, int n >
    inline void multiply ( const FieldMatrix< Field, m, n > &A,
                           const FieldVector< Field, n > &x,
                           FieldVector< Field, m > &y )
    {
      for( int i = 0; i < m; ++i )
      {
        Field &value = y[ i ];

        value = 0;
        for( int j = 0; j < n; ++j )
          value += A[ i ][ j ] * x[ j ];
      }
    }



    template< class Field, int m, int n, int p >
    inline void multiply ( const FieldMatrix< Field, m, n > &A,
                           const FieldMatrix< Field, n, p > &B,
                           FieldMatrix< Field, m, p > &C )
    {
      for( int i = 0; i < m; ++i )
      {
        for( int j = 0; j < p; ++j )
        {
          Field &value = C[ i ][ j ];
          
          value = 0;
          for( int k = 0; k < n; ++k )
            value += A[ i ][ k ] * B[ k ][ j ];
        }
      }
    }
    
  }
  
}

#endif
