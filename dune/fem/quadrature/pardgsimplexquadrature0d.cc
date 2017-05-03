#include <cassert>
#include <cmath>
#include <iostream>
#include <dune/fem/quadrature/pardgsimplexquadrature.hh>

namespace Dune {
namespace Fem {
namespace ParDGSimplexQuadrature
{

// return a quadrature formula that has at least degree minimum_degree
template<>
const Quadrature0d& Quadrature0d::quadrature(int minimum_degree)
{
  return quad0d;
}



static const int infinite = 1000000;
static const double quad0d_x[][1] = {{1.0}};
const Quadrature0d quad0d(1, infinite, quad0d_x);

} // ParDGSimplexQuadrature
} // Fem
} // Dune
