#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include <cassert>
#include "function.hpp"
#include "blas.hpp"


namespace pardg
{


template<int dim>
class Quadrature
{
public:
  Quadrature(int nop, int degree, const double x[][dim+1]);
  ~Quadrature();
  int number_of_points() const;
  const double* x(int i) const;
  double w(int i) const;
  const double* w() const;
  void integrate(Function &f, double *result) const;
  void check() const;

  static const Quadrature& quadrature(int minimum_degree);

private:
  typedef double Point[dim+1];

  const int nop;
  const int degree;
  Point *x_w;
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



} // namespace pardg



// class Quadrature<dim> inline implementation
template<int dim>
inline
pardg::Quadrature<dim>::Quadrature(int nop, int degree, 
				   const double x[][dim+1]) :
  nop(nop), degree(degree), x_w(new Point[nop])
{
  assert(x_w);
  
  for(int i=0; i<nop; i++){
    for(int l=0; l<=dim; l++) x_w[i][l] = x[i][l];
  }
}


template<int dim>
inline
pardg::Quadrature<dim>::~Quadrature()
{
  delete[] x_w;
}


template<int dim>
inline
int pardg::Quadrature<dim>::number_of_points() const
{
  return nop;
}


template<int dim>
inline
const double* pardg::Quadrature<dim>::x(int i) const
{
  assert(i>=0 && i<nop);
  return (dim == 0)? NULL : x_w[i];
}


template<int dim>
inline
double pardg::Quadrature<dim>::w(int i) const
{
  assert(i>=0 && i<nop);
  return x_w[i][dim];
}


template<int dim>
inline
const double* pardg::Quadrature<dim>::w() const
{
  return &x_w[0][dim];
}



template<int dim>
inline
void pardg::Quadrature<dim>::integrate(Function &f, double *result) const
{
  assert(dim == f.dim_of_argument());
  const int dim_of_value = f.dim_of_value();
  double *tmp = new double[dim_of_value];
  assert(tmp);

  dset(dim_of_value, 0.0, result, 1);
  for(int i=0; i<nop; i++){
    f(x(i), tmp);
    cblas_daxpy(dim_of_value, w(i), tmp, 1, result, 1);
  }

  delete[] tmp;
}



#endif
