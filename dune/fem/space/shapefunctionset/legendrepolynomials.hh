#ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_LEGENDREPOLYNOMIALS_HH
#define DUNE_FEM_SPACE_SHAPEFUNCTIONSET_LEGENDREPOLYNOMIALS_HH

// C++ inlcudes
#include <cassert>
#include <iostream>


namespace Dune
{

  namespace Fem
  {

    class LegendrePolynomials
    {
    public:
      static const int maxOrder = 11;

      static const double factor[ maxOrder ][ maxOrder ];
      static const double weight[ maxOrder ];

    public:
      static double evaluate ( const int num, const double x )
      {
        assert( 0 <= num && num < maxOrder );

        double phi = factor[ num ][ num ];
        for( int i = num-1; i >= 0; --i )
          phi = phi * x + factor[ num ][ i ];
        return weight[ num ] * phi;
      }

      static double jacobian ( const int num, const double x )
      {
        assert( 0 <= num && num < maxOrder );

        double phi = 0.;
        if( num >= 1 )
        {
          phi = factor[ num ][ num ] * num;
          for( int i = num-1; i >= 1; --i )
            phi = phi * x + factor[ num ][ i ] * i;
        }
        return weight[ num ] * phi;
      }

      static double hessian ( const int num, const double x )
      {
        assert( 0 <= num && num < maxOrder );

        double phi=0.;
        if( num >= 2 )
        {
          phi = factor[ num ][ num ] * num * ( num-1 );
          for( int i = num-1; i >= 2; --i )
            phi = phi * x + factor[ num ][ i ] * i * (i-1);
        }
        return weight[ num ] * phi;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_LEGENDREPOLYNOMIALS_HH
