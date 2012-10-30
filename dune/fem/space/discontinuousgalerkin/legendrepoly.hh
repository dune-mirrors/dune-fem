#ifndef DUNE_FEM_LEGENDREPOLY_HH
#define DUNE_FEM_LEGENDREPOLY_HH

#include <cassert>
#include <iostream>

namespace Dune
{

  namespace Fem 
  {  

    struct LegendrePoly
    {
      static double evaluate(const int num, const double x)
      {  
        assert(0<=num && num<maxPol);
        double phi = factor[ num ][ num ];
        for(int i=num-1; i>=0; --i )
        {
          phi = phi * x + factor[ num ][ i ];
        }
        return weight[ num ] * phi;
      }
      
      static double jacobian(const int num, const double x)
      { 
        assert(0<=num && num<maxPol);
        double phi=0.0;
        if (num>=1)
        {
          phi = factor[ num ][ num ] * num;
          for(int i=num-1 ; i>=1; --i )
          {
            phi = phi * x + factor[ num ][ i ] * i;
          }
        }
        return weight[ num ] * phi;
      }
         
      static double hessian(const int num, const double x)
      { 
        assert(0<=num && num<maxPol);
        double phi=0.0;    

        if (num >= 2)
        {
          phi = factor[ num ][ num ] * num * ( num-1 );
          for(int i = num-1; i>=2; --i )
          {
            phi=phi*x+factor[num][i]*i*(i-1);
          }
        }
        return weight[ num ] * phi;
      }
      
    protected:
      enum { maxPol = 11 };
      static const double factor[ maxPol ][ maxPol ];
      static const double weight[ maxPol ];
    };

  } // namespace Fem 

} // namespace Dune 

#endif // #ifndef DUNE_FEM_LEGENDREPOLY_HH
