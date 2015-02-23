#ifndef DUNE_FEM_QUAD_TEST_HH
#define DUNE_FEM_QUAD_TEST_HH

#include "../quadrature.hh"
//#include <dune/quadrature/fixedorder.hh>

namespace Dune {
  namespace Fem {

  class Quad_Test
  {
  public:
    virtual void run();

    template <class Quad, class Fixed>
    void fixedOrderComparisonExec(Quad& quad, Fixed& fixed);
    template <class Quad>
    void weightSummationExec(Quad& quad);
    template <class Quad>
    void integrationExec(Quad& quad);
    template <class Quad>
    void orderExec(Quad& quad, int order);
    template <class Quad>
    void sameGeometryExec(Quad& quad1, Quad& quad2);
    void indicesTest();

  private:
    struct Func {
      // Integral over reference line = -0.25
      double operator() (const FieldVector<double, 1>& x) {
        return 0.25 - x[0];
      }

      // Integral over reference triangle = -0.0333333333
      // Integral over reference quadrilateral = -0.27083333333
      double operator() (const FieldVector<double, 2>& x) {
        return
          (0.25 - x[0])*
          (0.1 + 0.7*x[0] + 0.8*x[1]);
      }

      // Integral over reference tetrahedron = -0.005513888888
      // Integral over reference hexahedron = -0.03854166666666
      // Integral over reference prism = -0.425
      // Integral over reference pyramid = ??
      double operator() (const FieldVector<double, 3>& x) {
        return
          (0.25 - x[0])*
          (0.1 + 0.7*x[0] + 0.8*x[1])*
          (x[0] - 2.0*x[1] + 0.9*x[2]);
      }
    };

  private:
  };

  } // end namespace Fem
} // end namespace Dune

#endif
