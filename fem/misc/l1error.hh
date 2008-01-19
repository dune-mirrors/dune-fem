/**************************************************************************
** Title: L1Error class
** $RCSfile$
** $Revision: 2630 $$Name$
** $Date: 2007-11-26 20:23:42 +0100 (Mon, 26 Nov 2007) $
** Copyright: GPL $Author: robertk $
** Description: L2 error class, which computes the error between a function
** and a discrete function. Extracted class from
** Roberts poisson-example
**
**************************************************************************/

#ifndef DUNE_L1ERROR_HH
#define DUNE_L1ERROR_HH

// where the quadratures are defined
#include <dune/fem/quadrature/cachequad.hh>

namespace Dune
{
  
/*======================================================================*/
/*!
 * \class L1Error
 * \brief The L1Error class provides methods for error computation
 *
 * The class calculates || u-u_h ||_L2
 *
 * Currently only error between a Function and a DiscreteFunction can be
 * computed or the error between two DiscreteFunctions.
 */
/*======================================================================*/
  
  template <class DiscreteFunctionType, int n=0>
  class L1Error
  {
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionType::RangeType RangeType;
    typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef typename GridType::template Codim<0>::Geometry EnGeometryType;
    typedef typename EntityType::ctype coordType;
    
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
      
    enum { dim = GridType::dimension};
    enum { spacePolOrd = DiscreteFunctionSpaceType :: polynomialOrder };
    enum { dimRange = RangeType :: dimension };

  public:
/*======================================================================*/
/*!
 * norm: computation of the error norm
 *
 * this method initiates the grid-walkthrough for error computation
 * The polynomial degree for numerical integration and the
 * Type of continuous function are passed as template arguments
 *
 * \param f the continuous function
 *
 * \param discFunc the discrete function
 *
 * \param time the time, at which the functions should
 * be evaluated
 *
 * \return the norm of the L2-Error as RangeType of DiscreteFunction
 */
/*======================================================================*/
    
    template <class FunctionType>
    RangeType norm (const FunctionType &f, DiscreteFunctionType &discFunc,
                    const double time)
    {
      return norm(f,discFunc,2*discFunc.space().order()+2,time);
    }

    template <class FunctionType>
    RangeType norm (const FunctionType &f, DiscreteFunctionType &discFunc,
                    int polOrd = (2 * spacePolOrd + 2),
                    double time = 0.0)
    {
      // get function space
      const DiscreteFunctionSpaceType & space = discFunc.space();
      
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
        // entity
        const EntityType& en = *it;
        
        // create quadrature for given geometry type
        CachingQuadrature <GridPartType , 0 > quad(en,polOrd);
        
        // get local function
        LocalFunctionType lf = discFunc.localFunction(en);

        // get geoemetry of entity
        const EnGeometryType& geo = en.geometry();
        
        // integrate
        const int quadNop = quad.nop();
        for(int qP = 0; qP < quadNop; ++qP)
        {
          const double det = quad.weight(qP) *
              geo.integrationElement(quad.point(qP));

          if (n==0) {
            f.evaluate(geo.global(quad.point(qP)),time,ret);
            lf.evaluate(quad[qP],phi);
            for(int k=0; k<dimRange; ++k)
            {
              error[k] += det * sqrt(SQR(ret[k] - phi[k]));
            }
          }
          else if (n==1)
          {
            f.evaluate(geo.global(quad.point(qP)),time,ret);
            lf.evaluate(quad[qP],phi);
            for(int k=0; k<dimRange; ++k) {
              error[k] += det * sqrt(SQR(ret[k] - phi[k]));
            }
          }
        }
      }
      
      for(int k=0; k<dimRange; ++k)
      {
        error[k] = comm.sum( error[k] );
      }
      
      return error;
    } // end of method

/*======================================================================*/
/*!
 * norm2: computation of difference between two discrete functions
 *
 * The last parameter dummy is in principle superfluous, as discrete
 * function do not know the time. But for analogy with the norm() method
 * the parameter is kept.
 *
 * \param f1 first discrete function
 *
 * \param f2 second discrete function
 *
 * \param dummy parameter, which can be empty (for continuous functions
 * this is used for passing the time)
 *
 * \return the l2-error
 */
/*======================================================================*/
  
  template <int polOrd>
  RangeFieldType norm2(const DiscreteFunctionType& f1,
               const DiscreteFunctionType& f2, double dummy = 0)
  {
    const DiscreteFunctionSpaceType & space = f1.space();
    
    const GridPartType & gridPart = space.gridPart();
    typedef typename GridPartType :: GridType :: Traits ::
        CollectiveCommunication
        CommunicatorType;
    
    const CommunicatorType & comm = gridPart.grid().comm();
    
    RangeFieldType ret=0;

    if( dimRange > 1 )
    {
      std::cout << "L1Error::norm2: only implemented for dimRange = 1! \n";
      abort();
    }
    
    // get function space
    const DiscreteFunctionSpaceType& dfsp = f1.space();
    
    RangeType lv1,lv2;
    
    // for product:
    //int quadOrd = (polOrd*dim)*(polOrd*dim);
    int quadOrd = polOrd;
                    
    // iterate over all elements defining the function
    IteratorType eit = dfsp.end();
    for (IteratorType it = dfsp.begin(); it!=eit; ++it)
    {
      const EntityType& en = *it;
      
      CachingQuadrature <GridPartType , 0 > quad(en,quadOrd);
      // get local functions on current element
      LocalFunctionType lf1 = f1.localFunction( en );
      LocalFunctionType lf2 = f2.localFunction( en );

      // get geoemetry of entity
      const EnGeometryType& geo = en.geometry();
      
      const int quadNop = quad.nop();
      for (int qp = 0; qp < quadNop; ++qp)
      {
        const double det =
            geo.integrationElement(quad.point(qp));

        // evaluate local functions
        lf1.evaluate(quad[qp], lv1);
        lf2.evaluate(quad[qp], lv2);
        // substract
        lv1 -= lv2;

        ret += det * quad.weight(qp) * sqrt(lv1 * lv1);
      } // end qp iteration
      
    } // end element iteration
    
    ret = comm.sum(ret);
    return ret;
  } // end method
  
}; // end of class L1Error

} // end namespace
#endif
