#ifndef DUNE_L2PROJECTION_HH
#define DUNE_L2PROJECTION_HH

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>

#include "dgl2projection.hh"

namespace Dune 
{

// deprecated class (use DGL2Projection instead)
struct L2ProjectionImpl : public DGL2ProjectionImpl
{
  //! project function onto discrete function space  
  template <class FunctionImp, class DiscreteFunctionImp>
  DUNE_DEPRECATED static void project(const FunctionImp& f, DiscreteFunctionImp& discFunc, int polOrd = -1)
  {
    DGL2ProjectionImpl :: project( f, discFunc, polOrd );
  }
};

/*======================================================================*/
/*! @ingroup L2ProjectionOperator
 *  \class L2Projection
 *  \brief The L2Projection class provides methods for projecting a function
 *         unto a given discrete function space. Note that this implementation
 *         assumes orthorgonal base functions!
 */
/*======================================================================*/
template <typename DFieldType, typename RFieldType,
          typename DType , typename RType>
class L2Projection : public Operator<DFieldType, RFieldType,DType , RType> {
 public:
  typedef DType DomainType;
  typedef RType  RangeType;
  typedef DFieldType DomainFieldType;
  typedef RFieldType RangeFieldType;

  //! Constructor taking degree for quadrature rule
  //! if no argument given a default value is chosen depending on the order
  //! in the discrete function space
  DUNE_VERSION_DEPRECATED(1,4,remove)
  L2Projection(int polOrd = -1) : polOrd_(polOrd) {}

  //! apply L2 projection
  virtual void operator() (const DomainType& f, RangeType& discFunc) const 
  {
    if( discFunc.space().continuous() )
      DUNE_THROW(NotImplemented,"L2-Projection not implemented for contiuous spaces!"); 
    else 
      DGL2ProjectionImpl::project(f,discFunc,polOrd_);
  }

protected:
  const int polOrd_;
};

  namespace Fem { 

    /*! @ingroup L2ProjectionOperator
     *  \class L2Projection
     *  \brief The L2Projection class provides methods for projecting a function
     *         unto a given discrete function space. Note that this implementation
     *         assumes orthorgonal base functions!
     */
    template < class DType, class RType>
    class L2Projection : public Operator<DType , RType> 
    {
    public:
      //! domain function type 
      typedef DType  DomainType;
      //! range function type 
      typedef RType  RangeType;

      /** \brief Constructor 
       *    \param  quadOrder      degree for quadrature rule (default = 2*space.order())
       *    \param  doCommunicate  apply communication for the result (default = true)
       */    
      explicit L2Projection( const int quadOrder = -1, 
                             const bool doCommunicate = true ) 
        : quadOrder_( quadOrder ), doCommunicate_( doCommunicate ) 
      {}

      /** \brief  calculates the L2 projection of a function onto the 
       *          discrete space discreteFunction belongs to.
       *  \param  function          function to be projected 
       *  \param  discreteFunction  discrete result of projection  */
      virtual void operator() ( const DomainType& function, RangeType& discreteFunction ) const 
      {
        if( discreteFunction.space().continuous() )
          DUNE_THROW(NotImplemented,"L2-Projection not implemented for contiuous spaces!"); 
        else 
          DGL2ProjectionImpl::project( function, discreteFunction, quadOrder_, doCommunicate_ );
      }

    protected:
      const int quadOrder_;          // order of quadrature  
      const bool doCommunicate_ ; // true if communication is applied for the result 
    };
  } // end namespace Fem 

}
#endif
