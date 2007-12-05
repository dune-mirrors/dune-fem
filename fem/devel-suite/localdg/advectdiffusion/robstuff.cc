#ifndef LOCALDG_STUFF_CC 
#define LOCALDG_STUFF_CC 

#include <dune/fem/quadrature/cachequad.hh>

namespace Dune {

// calculates || u-u_h ||_L2
template <class DiscreteFunctionType> 
class L2Error
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  
public:  
  template <class FunctionType> 
  RangeType norm (FunctionType &f, DiscreteFunctionType &discFunc,
      double time = 0.0)
  {
    const typename DiscreteFunctionType::FunctionSpaceType 
        & space = discFunc.space();  

    int polOrd = 2 * space.order() + 2;
  
    typedef typename FunctionSpaceType::GridType GridType;
    typedef typename FunctionSpaceType::GridPartType GridPartType;
    typedef typename FunctionSpaceType::IteratorType IteratorType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
   
    RangeType ret (0.0);
    RangeType phi (0.0);

    RangeType error(0.0);
    
    IteratorType it    = space.begin();
    IteratorType endit = space.end();

    // check whether grid is empty 
    assert( it != endit ); 
    
    CachingQuadrature <GridPartType , 0 > quad(*it,polOrd);
		
    enum { dimR = RangeType :: dimension };
    
    for(; it != endit ; ++it)
    {
      LocalFuncType lf = discFunc.localFunction(*it); 
      for(int qP = 0; qP < quad.nop(); qP++)
      {
        double det = quad.weight(qP) * (*it).geometry().integrationElement(quad.point(qP));
        f.evaluate(time,(*it).geometry().global(quad.point(qP)), ret);
        lf.evaluate( quad[ qP ], phi );
        for(int k=0; k<dimR; ++k)
          error[k] += det * SQR(ret[k] - phi[k]);
      }
    }
    
    for(int k=0; k<dimR; ++k)
      error[k] = sqrt(error[k]);

    return error;
  }
};

} // end namespace 

#endif

