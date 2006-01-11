#ifndef DUNE_PRISMPOINTS_HH
#define DUNE_PRISMPOINTS_HH

#include <dune/common/fvector.hh>

namespace Dune {
  
  //! Single point of reference for the quadrature points for prism.
  //! This class is a singleton, i.e. all points for all quadratures are
  //! created once.
  class PrismPoints
  {
  public:
    enum { numQuads = 3 };
    enum { MAXP=6 };
    enum { highest_order=2 };
    
    //! Access to the singleton object.
    inline static const PrismPoints& instance() {
      static PrismPoints prismPoints; 
      return prismPoints; 
    }

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
    
  private:
    //! initialize quadrature points on the interval for all orders
    PrismPoints() {
      int m = 0;
      O[m] = 0;
      N[m] = 0;
      
      // polynom degree 0  ???
      m = 1;
      G[m][0][0] = 0.0;
      G[m][0][1] = 0.0;
      G[m][0][2] = 0.0;
      
      G[m][1][0] = 1.0;
      G[m][1][1] = 0.0;
      G[m][1][2] = 0.0;
      
      G[m][2][0] = 0.0;
      G[m][2][1] = 1.0;
      G[m][2][2] = 0.0;
      
      G[m][3][0] = 0.0;
      G[m][3][1] = 0.0;
      G[m][3][2] = 1.0;
      
      G[m][4][0] = 1.0;
      G[m][4][1] = 0.0;
      G[m][4][2] = 1.0;
      
      G[m][5][0] = 0.0;
      G[m][5][1] = 0.1;
      G[m][5][2] = 1.0;

      W[m][0] = 0.08333333333333333;
      W[m][1] = 0.08333333333333333;
      W[m][2] = 0.08333333333333333;
      W[m][3] = 0.08333333333333333;
      W[m][4] = 0.08333333333333333;
      W[m][5] = 0.08333333333333333;
      
      O[m] = 0;// verify ????????
      N[m] = 6;
      
      // polynom degree 2  ???
      m = 2;
      G[m][0][0] =0.66666666666666666 ;
      G[m][0][1] =0.16666666666666666 ;
      G[m][0][2] =0.211324865405187 ;

      G[m][1][0] = 0.16666666666666666;
      G[m][1][1] =0.66666666666666666 ;
      G[m][1][2] = 0.211324865405187;
      
      G[m][2][0] = 0.16666666666666666;
      G[m][2][1] = 0.16666666666666666;
      G[m][2][2] = 0.211324865405187;
      
      G[m][3][0] = 0.66666666666666666;
      G[m][3][1] = 0.16666666666666666;
      G[m][3][2] = 0.788675134594813;
      
      G[m][4][0] = 0.16666666666666666;
      G[m][4][1] = 0.66666666666666666;
      G[m][4][2] = 0.788675134594813;
      
      G[m][5][0] = 0.16666666666666666;
      G[m][5][1] = 0.16666666666666666;
      G[m][5][2] = 0.788675134594813;
     
      W[m][0] = 0.08333333333333333;
      W[m][1] = 0.08333333333333333;
      W[m][2] = 0.08333333333333333;
      W[m][3] = 0.08333333333333333;
      W[m][4] = 0.08333333333333333;
      W[m][5] = 0.08333333333333333;

      O[m] = 2;// verify ????????
      N[m] = 6;
    }

  private:
    FieldVector<double, 3> G[numQuads][MAXP];//positions
    
    double W[numQuads][MAXP]; // weights associated with points       
    int O[numQuads];          // order of the rule
    int N[numQuads];          // number of points per quadrature rule
  };
}

#endif
