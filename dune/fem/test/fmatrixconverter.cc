#undef NDEBUG

#include <config.h>
#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/fem/misc/fmatrixconverter.hh>


using namespace Dune;
using namespace Fem;


const int n = 2;
const int m = 3;

typedef FieldVector< double, n*m > VectorType;
typedef FieldMatrix< double, n, m > MatrixType;

typedef FieldMatrixConverter< VectorType, MatrixType > ConverterType;

// main program
int main(int argc, char ** argv)
{
  try
  {
    VectorType vec( 4.0);
    MatrixType mat( 2.0 );

    ConverterType conv( vec );
    mat -= conv;


    return 0;
  }
  catch( Exception e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}
