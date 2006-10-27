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

#ifndef __DUNE_L2ERROR_HH__
#define __DUNE_L2ERROR_HH__

// where the quadratures are defined 
#include <dune/fem/quadrature/quadrature.hh>

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

template <class DiscreteFunctionType> 
class L2Error
{
  typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
  
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
 *   \param time the time, at which the functions should 
 *          be evaluated
 *
 *   \return the norm of the L2-Error
 */
/*======================================================================*/

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
   
    enum { dim = GridType :: dimension };
    Quadrature <typename FunctionSpaceType::RangeFieldType, dim> quad(
        it->geometry().type(), polOrd);
    
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
}; // end of class L2Error

} // end namespace 

#endif

