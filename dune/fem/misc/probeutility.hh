/**************************************************************************
**       Title: probeutility
**    $RCSfile: $
**   $Revision:$$Name:$
**       $Date:$
**   Copyright: GPL $Author:$
** Description: extract function values of a discrete function 
**************************************************************************/

#ifndef __PROBEUTILITY_HH__
#define __PROBEUTILITY_HH__

#include <cassert>
#include <cmath>

namespace Dune
{

/*======================================================================*/
/*!
 *  \class ProbeUtility
 *  \brief provides methods for probing of discretefunctions
 *
 *  Class collects routines for probing of discretefunctions. Currently
 *  point-probing and probing along a line are supported
 */
/*======================================================================*/

  class ProbeUtility
  {
  public:

/*======================================================================*/
/*!  @ingroup HelperClasses
 *   pointprobe: extract value at global point
 *
 *   \param func the discretefunction
 *
 *   \param point the point in global coordinates to be evaluated
 *
 *   \param value the return value
 */
/*======================================================================*/
    
    template <class DiscreteFunctionType>
    void pointprobe(DiscreteFunctionType& func, 
                    const typename DiscreteFunctionType::DomainType& point,
                    typename DiscreteFunctionType::RangeType& value)
          {
            typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType 
                DiscreteFunctionSpaceType;
            typedef typename DiscreteFunctionSpaceType::IteratorType 
                IteratorType;
            
            const DiscreteFunctionSpaceType& space = 
                func.space();
            
            // perform grid walkthrough and search point
            
            IteratorType endit = space.end();
            IteratorType it = space.begin();
            
            while( (it!=endit) && 
                   !(it->geometry().checkInside(it->geometry().local(point)))
                   )
                ++it;

            // return value in point            
            if (it==endit) // point not found
                value = NAN; // 
            else // iterator points to element containing point
            {
              typename DiscreteFunctionType::LocalFunctionType lf = 
                  func.localFunction(*it);
              lf.evaluate(it->geometry().local(point),value);
            }
          }


/*======================================================================*/
/*! 
 *   lineprobe: extract values along a line
 *
 *   the routine assumed a reasonable low number of points to be handled,
 *   as for all points, the local coordinates in all elements are computed
 *   the output array is assumed to be allocated sufficiently and must provide 
 *   a random access [] operator!!!
 *   So C-arrays and std::vectors, fielvectors, etc. can be used as 
 *   output types
 *
 *   \param func the discretefunction
 *
 *   \param start the start-point in global coordinates to be evaluated
 *
 *   \param end the end-point in global coordinates to be evaluated
 *
 *   \param npoints the number of points to be equally distributed 
 *          along the line 
 *
 *   \param values a pointer to an (allocated!) vector for the return values
 */
/*======================================================================*/
    
    template <class DiscreteFunctionType, class VectorType>
    void lineprobe(const DiscreteFunctionType& func, 
                   const typename DiscreteFunctionType::DomainType& start,
                   const typename DiscreteFunctionType::DomainType& end,
                   const int npoints,
                   VectorType& values)
          {
            // typedefs
            typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType 
                DiscreteFunctionSpaceType;
            typedef typename DiscreteFunctionType::DomainType DomainType;
            typedef typename DiscreteFunctionSpaceType::IteratorType 
                 IteratorType;
            
            assert(npoints>=2);
            
            // check, that values are allocated and initialize with nan
//            assert(values);
            for (int i=0;i<npoints;i++) values[i] = NAN;
            
            // field for marking, which points are finished
            bool pointfinished[npoints];
            for (int i=0;i<npoints;i++) pointfinished[i] = false;
            
            // field for storing the sample points
            DomainType points[npoints];
            DomainType step = end - start;
            step *= 1.0 /(npoints-1);
            DomainType p = start;
            for (int i=0;i<npoints;i++, p+=step)
                points[i] = p;
            
            const DiscreteFunctionSpaceType& space = 
                func.space();
            
            // perform grid walkthrough and search points            
            IteratorType endit = space.end();
            IteratorType it = space.begin();
            assert(it!=endit);
            
            // loop over all elements
            for (; it!=endit; ++it)
            {
              // loop over all points
              for (int i=0;i<npoints;i++)
                  // if not finished, check and compute values
                  if ( (!pointfinished[i]) &&
                       (it->geometry().checkInside(
                           it->geometry().local(points[i])))
                       )
                  {
                    typename DiscreteFunctionType::LocalFunctionType 
                        lf = func.localFunction(*it);
                    // lf.evaluate(it->geometry().local(points[i]),values[i]);
                    typename DiscreteFunctionType::RangeType val; 
                    lf.evaluate(it->geometry().local(points[i]),val);
                    values[i] = val;
                  }
              
            } // end of loop over grid
          } // end of lineprobe
    
  }; // end class ProbeUtility
  
} // end namespace Dune

#endif
