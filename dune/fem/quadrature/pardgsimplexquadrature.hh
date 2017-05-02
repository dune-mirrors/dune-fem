#ifndef DUNE_FEM_PARDGSIMPLEXQUADRATURE_HPP
#define DUNE_FEM_PARDGSIMPLEXQUADRATURE_HPP

#include <cassert>
#include <vector>

#include <dune/common/fvector.hh>

namespace Dune {
namespace Fem {
namespace ParDGSimplexQuadrature
{


template<int dim>
class Quadrature
{
public:
  enum { numCorners = dim+1 };
  typedef Dune::FieldVector< double, dim > CoordinateType;
  typedef Dune::FieldVector< double, dim+1 > Point;

  Quadrature(int nop, int degree, const double x[][dim+1])
    : nop(nop), degree(degree), x_w( nop )
  {
    for(int i=0; i<nop; i++)
    {
      for(int l=0; l<=dim; l++) x_w[i][l] = x[i][l];
    }
  }

  explicit Quadrature( const int order )
  {
    const Quadrature& quad = quadrature( order );
    nop    = quad.number_of_points();
    degree = quad.max_order();

    assert( order <= degree );

    x_w.resize( nop );
    for(int i=0; i<nop; i++)
    {
      for(int l=0; l<=dim; l++) x_w[i][l] = quad.x_w[i][l];
    }
  }

  CoordinateType point( const int i ) const
  {
    assert(i >= 0 && i < numPoints());
    CoordinateType result;
    for (size_t j = 0; j < dim; ++j)
    {
      result[j] = x(i)[j];
    }
    return result;
  }

  //! Access to the ith quadrature weight.
  double weight(const int i) const
  {
    assert(i >= 0 && i < numPoints());
    // scale with volume of reference element!
    return w(i);
  }

  int order() const { return degree; }

  int numPoints() const { return number_of_points(); }
  int max_order() const { return degree; }

  int number_of_points() const{ return nop; }

  const Point& x(int i) const {
    assert(i>=0 && i<nop);
    return x_w[i];
  }

  double w(int i) const
  {
    assert(i>=0 && i<nop);
    return x_w[i][dim];
  }

  void check() const;

  static const Quadrature& quadrature(int minimum_degree);

protected:

  int nop;
  int degree;
  std::vector< Point > x_w;
};


// dimensions of interest
typedef Quadrature<0> Quadrature0d;
typedef Quadrature<1> Quadrature1d;
typedef Quadrature<2> Quadrature2d;
typedef Quadrature<3> Quadrature3d;



// 0d quadrature rule
//
// there is only one with one point and infinite degree
extern const Quadrature0d quad0d;



// 1d quadrature rules for the interval (0,1)
//
// naming convention of quad1d_x: x denotes the order of the formula
//                                i.e. the maximum degree of polynimials
//                                that is integrated exact by the formula.
//
// exception quad1d_0: this is a dummy that has no points
extern const Quadrature1d quad1d_0, quad1d_1, quad1d_3, quad1d_5, quad1d_7,
  quad1d_9, quad1d_11, quad1d_13, quad1d_15, quad1d_17, quad1d_19, quad1d_21,
  quad1d_23, quad1d_25, quad1d_27, quad1d_29, quad1d_31, quad1d_33, quad1d_35,
  quad1d_37, quad1d_39;



// 2d quadrature rules for the reference triangle [(0,0), (1,0), (0,1)]
//
// naming convention of quad2d_x: x denotes the order of the formula
//                                i.e. the maximum degree of polynimials
//                                that is integrated exact by the formula.
//
// exception quad2d_0: this is a dummy that has no points
extern const Quadrature2d quad2d_0, quad2d_1, quad2d_2, quad2d_3, quad2d_4,
  quad2d_5, quad2d_6, quad2d_7, quad2d_8, quad2d_9, quad2d_10, quad2d_11,
  quad2d_13;




// 3d quadrature rules for the reference
// tetrahedron [(0,0,0), (1,0,0), (0,1,0), (0,0,1)]
//
// naming convention of quad3d_x: x denotes the order of the formula
//                                i.e. the maximum degree of polynimials
//                                that is integrated exact by the formula.
//
// exception quad3d_0: this is a dummy that has no points
extern const Quadrature3d quad3d_0, quad3d_1, quad3d_2, quad3d_3, quad3d_4,
  quad3d_5, quad3d_5b, quad3d_6, quad3d_7, quad3d_7b, quad3d_8, quad3d_9,
  quad3d_11;


} // ParDGSimplexQuadrature
} // Fem
} // Dune

#endif
