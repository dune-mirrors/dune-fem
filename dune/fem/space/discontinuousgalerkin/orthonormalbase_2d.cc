#ifndef DUNE_ORTHONORMALBASE_2D_CC
#define DUNE_ORTHONORMALBASE_2D_CC

#include <stdlib.h>
#include <stdio.h>
#include <cassert>

#include "orthonormalbase_2d.hh"

template class Dune :: Fem :: OrthonormalBase_2D< double, double >;
template class Dune :: Fem :: OrthonormalBase_2D< double, float  >;
template class Dune :: Fem :: OrthonormalBase_2D< float,  double >;
template class Dune :: Fem :: OrthonormalBase_2D< float,  float  >;
#endif
