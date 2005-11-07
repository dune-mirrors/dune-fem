#ifndef __DUNE_FEM_CC__
#define __DUNE_FEM_CC__

// assembling of the laplace operator using the 
// below defined getLocalMatrixMethod of the LaplaceOperator 

//#include <dune/fem/operator/laplace.hh>
#include "laplace.hh" 

// where the cg method is defined 
#include <dune/fem/inverseoperators.hh>

// where the quadratures are defined 
#include <dune/quadrature/fixedorder.hh>

namespace Dune 
{

// L2 projection of the rhs 
template <class DiscreteFunctionType> 
class L2Projection
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  
public:  
  template <int polOrd, class FunctionType> 
  void project (FunctionType &f, DiscreteFunctionType &discFunc)
  {
    const typename DiscreteFunctionType::FunctionSpaceType 
        & space = discFunc.getFunctionSpace();  
  
    discFunc.clear();
  
    typedef typename FunctionSpaceType::GridType GridType;
    typedef typename FunctionSpaceType::IteratorType IteratorType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
      

    typedef typename FunctionSpaceType::RangeType RangeType;
    RangeType ret (0.0);
    RangeType phi (0.0);

    IteratorType it    = space.begin();
    IteratorType endit = space.end();

    // check whether grid is empty 
    assert( it != endit ); 
    
    FixedOrderQuad <typename FunctionSpaceType::RangeFieldType,
              typename FunctionSpaceType::DomainType , polOrd > quad ( *it );
              
    //LocalFuncType lf = discFunc.newLocalFunction(); 
    
    for( ; it != endit ; ++it)
    {
      //discFunc.localFunction( *it , lf ); 
      LocalFuncType lf = discFunc.localFunction( *it ); 

      const typename FunctionSpaceType::BaseFunctionSetType & set = 
            space.getBaseFunctionSet(*it);

      for(int i=0; i<lf.numDofs(); i++)
      {
        for(int qP = 0; qP < quad.nop(); qP++)
        {
          double det = (*it).geometry().integrationElement(quad.point(qP));
          f.evaluate((*it).geometry().global( quad.point(qP) ), ret);
          set.eval(i,quad,qP,phi);
          lf[i] += det * quad.weight(qP) * (ret * phi);
        }
      }
    }
  }

};

// used for calculation of the initial values 
template <class DiscreteFunctionType> 
class LagrangeInterpolation
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  
public:  
  template <class FunctionType> 
  void interpol (FunctionType &f, DiscreteFunctionType &discFunc)
  {
    const typename DiscreteFunctionType::FunctionSpaceType
        & space = discFunc.getFunctionSpace();  
  
    typedef typename FunctionSpaceType::GridType GridType;
    typedef typename FunctionSpaceType:: IteratorType IteratorType; 
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
      
    typedef typename FunctionSpaceType::RangeType RangeType;
    RangeType ret (0.0);
    //LocalFuncType lf = discFunc.newLocalFunction(); 

    IteratorType endit = space.end(); 
    for(IteratorType it = space.begin(); it != endit ; ++it)
    {
      //discFunc.localFunction(*it,lf); 
      LocalFuncType lf = discFunc.localFunction(*it); 
      int numDof = lf.numDofs ();  
      for(int i=0; i<numDof; i++)
      {
        f.evaluate(it->geometry()[i],ret);
        lf[i] = ret[0];
      }
    }
  }
};

// calculates || u-u_h ||_L2
template <class DiscreteFunctionType> 
class L2Error
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  
public:  
  template <int polOrd, class FunctionType> 
  double norm (FunctionType &f, DiscreteFunctionType &discFunc,
      double time)
  {
    const typename DiscreteFunctionType::FunctionSpaceType 
        & space = discFunc.getFunctionSpace();  
  
    typedef typename FunctionSpaceType::GridType GridType;
    typedef typename FunctionSpaceType::IteratorType IteratorType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
   
    typedef typename FunctionSpaceType::RangeType RangeType;
    
    RangeType ret (0.0);
    RangeType phi (0.0);

    double sum = 0.0;
    //LocalFuncType lf = discFunc.newLocalFunction(); 
    
    IteratorType it    = space.begin();
    IteratorType endit = space.end();

    // check whether grid is empty 
    assert( it != endit ); 
    
    FixedOrderQuad < typename FunctionSpaceType::RangeFieldType,
               typename FunctionSpaceType::DomainType , polOrd > quad ( *it );
    
    for(; it != endit ; ++it)
    {
      LocalFuncType lf = discFunc.localFunction(*it); 
      for(int qP = 0; qP < quad.nop(); qP++)
      {
        double det = (*it).geometry().integrationElement(quad.point(qP));
        f.evaluate((*it).geometry().global(quad.point(qP)),time, ret);
        lf.evaluate((*it),quad,qP,phi);
        sum += det * quad.weight(qP) * SQR(ret[0] - phi[0]);
      }
    }
    return sqrt(sum);
  }
};

} // end namespace 

#endif

