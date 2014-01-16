#ifndef DUNE_ORTHONORMALBASE_3D_CC
#define DUNE_ORTHONORMALBASE_3D_CC

#include <stdlib.h>
#include <stdio.h>
#include <cassert>

#include "orthonormalbase_3d.hh"

template class Dune :: Fem :: OrthonormalBase_3D< double, double >;
template class Dune :: Fem :: OrthonormalBase_3D< double, float  >;
template class Dune :: Fem :: OrthonormalBase_3D< float,  double >;
template class Dune :: Fem :: OrthonormalBase_3D< float,  float  >;
#endif
