#ifndef DUNE_FEM_PYRAMIDPOINTS_HH
#define DUNE_FEM_PYRAMIDPOINTS_HH

#include <dune/common/fvector.hh>
#include <dune/common/visibility.hh>

namespace Dune
{

  namespace Fem
  {

    //! Single point of reference for the quadrature points for pyramids.
    //! This class is a singleton, i.e. all points for all quadratures are
    //! created once.
    class PyramidPoints
    {
    public:
      enum { numQuads = 2 };
      enum { MAXP = 8 };
      enum { highest_order = 2 };

      //! Access to the ith point of quadrature rule m.
      const FieldVector<double, 3>& point(int m, int i) const
      {
        return G[m][i];
      }

      //! Access to the ith weight of quadrature rule m.
      double weight (int m, int i) const
      {
        return W[m][i];
      }

      //! Actual order of quadrature rule m.
      int order (int m) const
      {
        return O[m];
      }

      //! Number of points in the quadrature rule m.
      int numPoints(int m) const
      {
        return N[m];
      }

      //! initialize quadrature points on the interval for all orders
      PyramidPoints()
      {
        int m = 0;
        O[m] = 0;
        N[m] = 0;

        // polynom degree 2  ???
        m = 1;
        G[m][0][0] =0.58541020;
        G[m][0][1] =0.72819660;
        G[m][0][2] =0.13819660;

        G[m][1][0] =0.13819660;
        G[m][1][1] =0.72819660;
        G[m][1][2] =0.13819660;

        G[m][2][0] =0.13819660;
        G[m][2][1] =0.27630920;
        G[m][2][2] =0.58541020;

        G[m][3][0] =0.13819660;
        G[m][3][1] =0.27630920;
        G[m][3][2] =0.13819660;

        G[m][4][0] =0.72819660;
        G[m][4][1] =0.13819660;
        G[m][4][2] =0.13819660;

        G[m][5][0] =0.72819660;
        G[m][5][1] =0.58541020;
        G[m][5][2] =0.13819660;

        G[m][6][0] =0.27630920;
        G[m][6][1] =0.13819660;
        G[m][6][2] =0.58541020;

        G[m][7][0] =0.27630920;
        G[m][7][1] =0.13819660;
        G[m][7][2] =0.13819660;

        W[m][0] = 0.041666667;
        W[m][1] = 0.041666667;
        W[m][2] = 0.041666667;
        W[m][3] = 0.041666667;
        W[m][4] = 0.041666667;
        W[m][5] = 0.041666667;
        W[m][6] = 0.041666667;
        W[m][7] = 0.041666667;

        O[m] = 2;// verify ????????
        N[m] = 8;
      }

    private:
      FieldVector<double, 3> G[numQuads][MAXP];
      double W[numQuads][MAXP]; // weights associated with points
      int O[numQuads];          // order of the rule
      int N[numQuads];          // number of points in quadrature rule
    };

  } // namespace Fem

}  // namespace Dune

#endif // #ifndef DUNE_FEM_PYRAMIDPOINTS_HH
