/**************************************************************************
**       Title: L2Error class
**    $RCSfile$
**   $Revision$$Name$
**       $Date$
**   Copyright: GPL $Author$
** Description: L2 error class, which computes the error between a function
**              and a discrete function. Extracted class from
**              Roberts poisson-example 
**
**************************************************************************/

#ifndef DUNE_L2ERROR_HH
#define DUNE_L2ERROR_HH

// where the quadratures are defined 
#include <dune/fem/quadrature/cachequad.hh>
#include <dune/grid/utility/twistutility.hh>

namespace Dune 
{
  
/*======================================================================*/
/*!
 *  \class L2Error
 *  \brief The L2Error class provides methods for error computation
 *
 *  The class calculates || u-u_h ||_L2
 *
 *  Currently only error between a Function and a DiscreteFunction can be
 *  computed or the error between two DiscreteFunctions. 
 */
/*======================================================================*/
  
  template <class DiscreteFunctionType,int n=0> 
  class L2Error
  {
    typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
    typedef typename DiscreteFunctionType::FunctionSpaceType DiscFuncSpaceType;
    typedef typename DiscreteFunctionType::RangeType         RangeType;
    typedef typename DiscFuncSpaceType::IteratorType         IteratorType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
    typedef typename DiscFuncSpaceType::GridType             GridType;
    typedef typename DiscFuncSpaceType::GridPartType         GridPartType;
    typedef typename GridType::template Codim<0>::Entity     EntityType; 
    typedef typename EntityType::ctype                       coordType; 
    
    enum {dimrange = DiscFuncSpaceType::DimRange};
    enum {dim = GridType::dimension};
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType 
    DiscreteFunctionSpaceType;
    enum { spacePolOrd = DiscreteFunctionSpaceType :: polynomialOrder }; 
    enum { dimR = RangeType :: dimension };

  public:  
/*======================================================================*/
/*! 
 *   norm: computation of the error norm
 *
 *   this method initiates the grid-walkthrough for error computation
 *   The polynomial degree for numerical integration and the 
 *   Type of continuous function are passed as template arguments
 *
 *   \param f the continuous function
 *
 *   \param discFunc the discrete function
 *
 *   \param polOrd for Quadrature to use, default is 2 * spaceOrd + 2 
 *   \param time the time, at which the functions should 
 *          be evaluated
 *
 *   \return the norm of the L2-Error as RangeType of DiscreteFunction 
 */
/*======================================================================*/
    
    template <class FunctionType> 
    RangeType norm (const FunctionType &f, DiscreteFunctionType &discFunc,
                    int polOrd = (2 * spacePolOrd + 2), 
                    double time = 0.0)
          {
            const DiscreteFunctionSpaceType & space = discFunc.space();  
            
            typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
            typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
            typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;
            
            const GridPartType & gridPart = space.gridPart();
            typedef typename GridPartType :: GridType :: Traits :: 
                CollectiveCommunication
                CommunicatorType; 
            
            const CommunicatorType & comm = gridPart.grid().comm();
            
            RangeType ret (0.0);
            RangeType phi (0.0);
            
            RangeType error(0.0);
            
            IteratorType endit = space.end();
            for(IteratorType it = space.begin(); it != endit ; ++it)
            {
              // create quadrature for given geometry type 
              CachingQuadrature <GridPartType , 0 > quad(*it,polOrd); 
              // get local function 
              LocalFuncType lf = discFunc.localFunction(*it); 
              
              // integrate 
              const int quadNop = quad.nop();
              for(int qP = 0; qP < quadNop; ++qP)
              {
                double det = quad.weight(qP) * 
                    (*it).geometry().integrationElement(quad.point(qP));
                if (n==0) {
                  f.evaluate((*it).geometry().global(quad.point(qP)),time,ret);
                  lf.evaluate(quad,qP,phi);
                  for(int k=0; k<dimR; ++k) {
                    error[k] += det * SQR(ret[k] - phi[k]);
                  }
                } else if (n==1) {
                  
                  f.evaluate((*it).geometry().global(quad.point(qP)),time,ret);
                  lf.evaluate(quad,qP,phi);
                  for(int k=0; k<dimR; ++k) {
                    error[k] += det * SQR(ret[k] - phi[k]);
                  }
                }
              }
            }
            
            //return sqrt(sum);
            
            for(int k=0; k<dimR; ++k)
            {
              error[k] = comm.sum( error[k] );
              error[k] = sqrt(error[k]);
            }
            
            return error;
          } // end of method

/*======================================================================*/
/*! 
 *   norm2: computation of difference between two discrete functions
 *
 *   The last parameter dummy is in principle superfluous, as discrete 
 *   function do not know the time. But for analogy with the norm() method
 *   the parameter is kept.
 *
 *   \param f1 first discrete function 
 *    
 *   \param f2 second discrete function 
 *
 *   \param dummy parameter, which can be empty (for continuous functions 
 *          this is used for passing the time) 
 *
 *   \return the l2-error
 */
/*======================================================================*/
  
  template <int polOrd> 
  double norm2(const DiscreteFunctionType& f1,
               const DiscreteFunctionType& f2, double dummy = 0)
  {
    const DiscreteFunctionSpaceType & space = f1.space();  
    
    const GridPartType & gridPart = space.gridPart();
    typedef typename GridPartType :: GridType :: Traits :: 
        CollectiveCommunication
        CommunicatorType; 
    
    const CommunicatorType & comm = gridPart.grid().comm();
    
    double ret=0;  
    
    assert( dimrange == 1); // currently only scalar functions supported
    
    const DiscFuncSpaceType& dfsp = f1.getFunctionSpace();  
    
    IteratorType it = dfsp.begin();  
    IteratorType eit = dfsp.end();  
    RangeType lv1,lv2;
    
    // get quadrature rule for first entity and hope that the geometry 
    // does not change. Get rule which is exact for product of df1,df2, 
    // assuming them to be polynomial
    
    int quadOrd = (polOrd*dim)*(polOrd*dim);
                    
    // iterate over all elements defining the function
    for (IteratorType it = dfsp.begin(); it!=eit; ++it)
    {
      CachingQuadrature <GridPartType , 0 > quad(*it,quadOrd); 
      // get local functions on current element
      LocalFunctionType lf1 = f1.localFunction( *it ); 
      LocalFunctionType lf2 = f2.localFunction( *it );
      for (int qp = 0; qp < quad.nop(); qp++)
      {
        const double det = 
            it->geometry().integrationElement(quad.point(qp));
        
      // the following demonstrates a very strange effect:
        lf1.evaluate(quad.point(qp), lv1);
        lf2.evaluate(quad.point(qp), lv2);
        ret += det * quad.weight( qp) * (lv1 - lv2) * (lv1 - lv2);
      } // end qp iteration
      
    } // end element iteration
    
    ret = comm.sum(ret);
    return sqrt(ret);
  }; // end method
  
}; // end of class L2Error

} // end namespace 
#endif
