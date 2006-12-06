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
 *  computed
 */
/*======================================================================*/

  template <class DiscreteFunctionType,int n=0> 
class L2Error
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

  enum { spacePolOrd = DiscreteFunctionSpaceType :: polynomialOrder };
  
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
    typedef typename GridPartType :: GridType :: Traits :: CollectiveCommunication
      CommunicatorType; 

    const CommunicatorType & comm = gridPart.grid().comm();
   
    RangeType ret (0.0);
    RangeType phi (0.0);
    
    RangeType error(0.0);

    enum { dimR = RangeType :: dimension };
    
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
        double det = quad.weight(qP) * (*it).geometry().integrationElement(quad.point(qP));
	if (n==0) {
	  f.evaluate((*it).geometry().global(quad.point(qP)),time,ret);
	  lf.evaluate(quad,qP,phi);
	  for(int k=0; k<dimR; ++k) {
	    error[k] += det * SQR(ret[k] - phi[k]);
	  }
	} else if (n==1) {
	  
	  f.evaluate((*it).geometry().global(quad.point(qP)),time,ret);
	  lf.evaluate((*it),quad,qP,phi);
	  for(int k=0; k<dimR; ++k) {
	    error[k] += det * SQR(ret[k] - phi[k]);
	  }
	}
      }
    }
    
    //comm.sum( &error[0] , dimR );

    for(int k=0; k<dimR; ++k)
    {
      error[k] = comm.sum( error[k] );
      error[k] = sqrt(error[k]);
    }

    return error;
  }
}; // end of class L2Error

} // end namespace 
#endif
