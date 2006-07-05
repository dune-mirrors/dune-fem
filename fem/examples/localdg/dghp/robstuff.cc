#ifndef LOCALDG_STUFF_CC 
#define LOCALDG_STUFF_CC 

#include <../../../quadrature/cachequad.hh>
#include <../../../quadrature/gausspoints.hh>

namespace Dune {

// calculates || u-u_h ||_L1
template <class DiscreteFunctionType> 
class L1Error
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  
public:  
  template <class FunctionType> 
  RangeType norm (FunctionType &f, DiscreteFunctionType &discFunc,
      double time = 0.0)
  {
    const typename DiscreteFunctionType::FunctionSpaceType 
        & space = discFunc.getFunctionSpace();  

    int polOrd = 2 * space.polynomOrder() + 2;
  
    typedef typename FunctionSpaceType::GridType GridType;
    typedef typename FunctionSpaceType::IteratorType IteratorType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
   
    RangeType ret (0.0);
    RangeType phi (0.0);

    RangeType error(0.0);
    
    IteratorType it    = space.begin();
    IteratorType endit = space.end();

    // check whether grid is empty 
    assert( it != endit ); 
    
    CachingQuadrature <GridType , 0 > quad(*it,polOrd);
		
    enum { dimR = RangeType :: dimension };
    
    for(; it != endit ; ++it)
    {
      LocalFuncType lf = discFunc.localFunction(*it); 
      for(int qP = 0; qP < quad.nop(); qP++)
      {
        double det = 
	  quad.weight(qP) * (*it).geometry().integrationElement(quad.point(qP));
        f.evaluate(time,(*it).geometry().global(quad.point(qP)), ret);
        lf.evaluate((*it),quad,qP,phi);
        for(int k=0; k<dimR; ++k)
          error[k] += det * fabs(ret[k] - phi[k]);
      }
    }
    
    //for(int k=0; k<dimR; ++k)
    //  error[k] = sqrt(error[k]);

    return error;
  }
};
// calculates || u-u_h ||_L1
template <class DiscreteFunctionType> 
class L1L1Error
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  typedef typename DiscreteFunctionType::RangeType RangeType;
  
public:  
  template <class FunctionType> 
  RangeType norm (FunctionType &f, DiscreteFunctionType &discFunc,
		  double time = 0.0,double dt = 1.0)
  {
    const typename DiscreteFunctionType::FunctionSpaceType 
        & space = discFunc.getFunctionSpace();  

    int polOrd = 2 * space.polynomOrder() + 2;
  
    typedef typename FunctionSpaceType::GridType GridType;
    typedef typename FunctionSpaceType::IteratorType IteratorType;
   
    RangeType ret (0.0);
    RangeType phi (0.0);

    RangeType error(0.0);
    
    IteratorType it    = space.begin();
    IteratorType endit = space.end();

    // check whether grid is empty 
    assert( it != endit ); 
    
    CachingQuadrature <GridType , 0 > quad(*it,polOrd);
		
    enum { dimR = RangeType :: dimension };
    
    const int nopTime = 3; 
    const GaussPts& timeQuad = GaussPts::instance();

    for(; it != endit ; ++it) {
      discFunc.setEntity(*it);
      for(int qP = 0; qP < quad.nop(); qP++) {
        double det = 
	  dt*quad.weight(qP) * 
	  (*it).geometry().integrationElement(quad.point(qP));
	for (int qT=0;qT<nopTime;++qT) {
	  double s = timeQuad.point(nopTime,qT);
	  double ldet = det*timeQuad.weight(nopTime,qT);
	  f.evaluate(time+s*dt,(*it).geometry().global(quad.point(qP)), ret);
	  // lf.evaluate((*it),quad,qP,phi);
	  phi = discFunc.uval(*it,quad,qP,s,2);
	  for(int k=0; k<dimR; ++k)
	    error[k] += ldet * fabs(ret[k] - phi[k]);
	}
      }
    }
    
    //for(int k=0; k<dimR; ++k)
    //  error[k] = sqrt(error[k]);

    return error;
  }
};

} // end namespace 

#endif

