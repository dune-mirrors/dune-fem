#ifndef DUNE_ORTHONORMALBASE_1D_CC
#define DUNE_ORTHONORMALBASE_1D_CC

#include <stdlib.h>
#include <stdio.h>
#include <cassert>

#include "orthonormalbase_1d.hh"

template class Dune :: Fem :: OrthonormalBase_1D< double, double >;
template class Dune :: Fem :: OrthonormalBase_1D< double, float  >;
template class Dune :: Fem :: OrthonormalBase_1D< float,  double >;
template class Dune :: Fem :: OrthonormalBase_1D< float,  float  >;
#endif
